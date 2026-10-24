#!/usr/bin/env python3
"""Which tests a change needs, for CI.

    ci_select.py CHANGED_FILES BUILD_DIR OUT_TESTS

CHANGED_FILES lists the paths the push or pull request changed, one per line, or holds the word
ALL. BUILD_DIR is the configured (not yet built) build: the tests' commands come from its
CTestTestfile.cmake, since `ctest --show-only` leaves out a command whose program is not built yet,
and their fixtures from `ctest --show-only=json-v1`. Writes the selected test
names, one per line, to OUT_TESTS, and prints to stdout, for $GITHUB_OUTPUT:

    mode=all|subset|none
    targets=<the CMake targets those tests and their fixtures' setups run>

A change to the library, the build, the CI or anything this script does not know runs everything.
A change to one binary's own source (src/niceops.cpp, tools/<name>.cpp, tests/check_*.cpp, the
Catch2 sources) runs the tests that call that binary, the tests reading what those write (fixtures,
followed transitively), and builds only the targets they need. A golden file runs the tests that
name it. Documentation runs nothing (the workflow's paths-ignore usually stops it earlier).
"""
import json
import re
import subprocess
import sys

IGNORED = [
    r'^docs/', r'\.md$', r'^reports/', r'^research_notes/', r'^LICENSE$', r'^[^/]*batch\.sh$',
    r'^tools/[^/]*\.py$', r'^tests/[^/]*\.py$', r'^\.claude/', r'^\.gitignore$',
]


def executables(cmake_lists):
    """{source path: target} for every add_executable whose first source is its own file."""
    found = {}
    for m in re.finditer(r'add_executable\(\s*([^\s)]+)\s+([^\s)]+)', open(cmake_lists).read()):
        target, source = m.group(1), m.group(2)
        if source.endswith('.cpp'):
            found[source] = target
    return found


def binaries_for(path, own_sources, unit_target):
    """The targets a change to `path` touches, None for everything, or [] for nothing."""
    if any(re.search(p, path) for p in IGNORED):
        return []
    if path in own_sources:
        return [own_sources[path]]
    if re.match(r'^tests/test_[^/]*\.cpp$', path) or re.match(r'^tests/[^/]*\.hpp$', path):
        return [unit_target]
    return None


def load_tests(build):
    """The tests, each with its expanded command (from CTestTestfile.cmake) and properties."""
    listing = subprocess.run(['ctest', '--test-dir', build, '--show-only=json-v1'],
                             capture_output=True, text=True, check=True).stdout
    tests = json.loads(listing)['tests']
    commands = {}
    for line in open(f'{build}/CTestTestfile.cmake'):
        m = re.match(r'add_test\((\S+) (.*)\)\s*$', line)
        if m:
            commands[m.group(1).strip('[]=')] = re.findall(r'"((?:[^"\\]|\\.)*)"', m.group(2))
    for t in tests:
        t['command'] = commands.get(t['name'], t.get('command', []))
        if not t['command']:
            sys.exit(f'ci_select: no command found for {t["name"]}')
    return tests


def main():
    changed_path, build, out_path = sys.argv[1:4]
    changed = [line.strip() for line in open(changed_path) if line.strip()]
    tests = load_tests(build)

    def props(t, key):
        for p in t.get('properties', []):
            if p['name'] == key:
                return p['value'] if isinstance(p['value'], list) else [p['value']]
        return []

    def write(names, mode, targets):
        with open(out_path, 'w') as f:
            f.write(''.join(n + '\n' for n in names))
        print(f'mode={mode}')
        print('targets=' + ' '.join(sorted(targets)))

    if changed == ['ALL']:
        write([t['name'] for t in tests], 'all', [])
        return

    own = executables('CMakeLists.txt')
    unit = own.get('tests/test_build.cpp', 'niceops_unit_tests')
    binaries, goldens = set(), set()
    for path in changed:
        if path.startswith('golden/'):
            goldens.add(path)
            continue
        hit = binaries_for(path, own, unit)
        if hit is None:
            write([t['name'] for t in tests], 'all', [])
            sys.stderr.write(f'{path} touches everything\n')
            return
        binaries.update(hit)
    if not binaries and not goldens:
        write([], 'none', [])
        return

    all_targets = set(own.values()) | {'niceops_unit_tests', 'niceops_plotting_stub'}

    def uses(t):
        """The CMake targets a test's command runs or reads."""
        found = set()
        for arg in t['command']:
            for target in all_targets:
                if re.search(r'(^|[/\s;="])' + re.escape(target) + r'($|[\s;"])', arg) \
                        or f'/CMakeFiles/{target}.dir/' in arg:
                    found.add(target)
        return found

    # A hidden Catch2 case ([.tag]) writes input files for CTests; it is the changed binary's test
    # only where a changed source holds its tag.  Otherwise the unit tests are the shards.
    changed_text = ''.join(open(p).read() for p in changed if p.startswith('tests/') and p.endswith(('.cpp', '.hpp')))

    def calls(t):
        hit = uses(t) & binaries
        if hit == {unit}:
            tags = [a for a in t['command'] if a.startswith('[.')]
            return not tags or any(tag in changed_text for tag in tags)
        return bool(hit)

    by_name = {t['name']: t for t in tests}
    selected = {t['name'] for t in tests
                if calls(t) or any(g in arg for g in goldens for arg in t['command'])}
    # What reads the selected runs' output: every test requiring a fixture they set up.
    while True:
        made = {f for n in selected for f in props(by_name[n], 'FIXTURES_SETUP')}
        more = {t['name'] for t in tests if set(props(t, 'FIXTURES_REQUIRED')) & made} - selected
        if not more:
            break
        selected |= more
    # What ctest will add: the setups of every fixture the selection requires, transitively.
    run = set(selected)
    while True:
        needed = {f for n in run for f in props(by_name[n], 'FIXTURES_REQUIRED')}
        more = {t['name'] for t in tests if set(props(t, 'FIXTURES_SETUP')) & needed} - run
        if not more:
            break
        run |= more
    targets = set().union(*(uses(by_name[n]) for n in run)) if run else set()
    write(sorted(selected), 'subset' if selected else 'none', targets)


if __name__ == '__main__':
    main()

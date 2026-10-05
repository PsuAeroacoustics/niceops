import re, sys
L = open('/Users/evg5332/.claude/jobs/a090ec0c/tmp/aamp.txt').read().splitlines()
for phrase in sys.argv[1:]:
    page = None; hits = []
    for line in L:
        m = re.match(r'=====PDFPAGE (\d+)', line)
        if m: page = int(m.group(1)); continue
        if phrase.lower() in line.lower(): hits.append(page)
    print(f"{phrase!r}: pdf {sorted(set(hits))} -> printed {[h - 16 for h in sorted(set(hits))]}")

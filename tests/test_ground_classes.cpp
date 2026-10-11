// The ground classes reader (--ground_classes): its columns, inheritance from
// the command line's ground, --units, and every refusal, on files written to a
// scratch folder.  What a run does with the classes is checked by the
// ground_class* CTests.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>

#include "ground_classes.hpp"

using namespace niceops;
using Catch::Matchers::ContainsSubstring;

namespace {

// A scratch folder holding one classes file, removed with everything in it.
struct Folder
{
    std::filesystem::path root;
    Folder()
    {
        char pattern[] = "/tmp/niceops_ground_classes_XXXXXX";
        REQUIRE(::mkdtemp(pattern) != nullptr);
        root = pattern;
    }
    ~Folder()
    {
        std::error_code ec;
        std::filesystem::remove_all(root, ec);
    }
    std::string write(const std::string& contents) const
    {
        const auto path = root / "classes.csv";
        std::ofstream(path) << contents;
        return path.string();
    }
};

// A command-line ground unlike every default, so inheritance shows.
GroundModel run_ground()
{
    GroundModel ground;
    ground.flow_resistance = 200.0;
    ground.impedance_model = ImpedanceModel::variable_porosity;
    ground.porosity_rate = 10.0;
    ground.layer_depth = 0.1;
    ground.roughness = 0.02;
    ground.receiver_height = 3.0;
    ground.turbulence = 5e-4;
    ground.turbulence_model = TurbulenceModel::harmonoise;
    ground.gamma_t = 3e-5;
    return ground;
}

} // namespace

TEST_CASE("a ground classes file reads its columns by their normalized names", "[ground_classes]")
{
    Folder folder;
    const auto path = folder.write("# surfaces\n"
                                   "\n"
                                   "Class, Name, Impedance_Model, Flow_Resistance (kPa s/m^2), Porosity_Rate (1/m), "
                                   "Layer_Depth (m), Roughness (m)\r\n"
                                   "0, grass, variable_porosity, 200, 10, , 0.02\n"
                                   "# a row left out\n"
                                   "21, pad, delany_bazley_layer, 10000, , 0.05, 0\n");
    const GroundClassFile file = read_ground_classes(path);
    REQUIRE(file.rows.size() == 2);
    REQUIRE(file.layer_depth_unit == "m");
    REQUIRE(file.roughness_unit == "m");
    REQUIRE(file.columns.count("flowresistance") == 1);
    const GroundClassRow& grass = file.rows[0];
    REQUIRE(grass.line == 4);
    REQUIRE(grass.id == 0);
    REQUIRE(grass.name == "grass");
    REQUIRE(grass.impedance_model == ImpedanceModel::variable_porosity);
    REQUIRE(grass.flow_resistance == 200.0);
    REQUIRE(grass.porosity_rate == 10.0);
    REQUIRE_FALSE(grass.layer_depth);
    REQUIRE(grass.roughness == 0.02);
    const GroundClassRow& pad = file.rows[1];
    REQUIRE(pad.line == 6);
    REQUIRE(pad.id == 21);
    REQUIRE(pad.impedance_model == ImpedanceModel::delany_bazley_layer);
    REQUIRE_FALSE(pad.porosity_rate);
    REQUIRE(pad.layer_depth == 0.05);
    REQUIRE(pad.roughness == 0.0);
}

TEST_CASE("a ground class takes what its row leaves out from the command line", "[ground_classes]")
{
    Folder folder;
    const GroundModel ground = run_ground();
    SECTION("a row of only an id is the command line's ground, field for field")
    {
        const auto classes = ground_classes_from(read_ground_classes(folder.write("class\n7\n")), ground, meters);
        REQUIRE(classes.size() == 1);
        const GroundModel& model = classes[0].model;
        REQUIRE(classes[0].id == 7);
        REQUIRE(classes[0].name == "7");
        REQUIRE(model.flow_resistance == ground.flow_resistance);
        REQUIRE(model.impedance_model == ground.impedance_model);
        REQUIRE(model.porosity_rate == ground.porosity_rate);
        REQUIRE(model.layer_depth == ground.layer_depth);
        REQUIRE(model.roughness == ground.roughness);
        REQUIRE(model.receiver_height == ground.receiver_height);
        REQUIRE(model.turbulence == ground.turbulence);
        REQUIRE(model.turbulence_model == ground.turbulence_model);
        REQUIRE(model.gamma_t == ground.gamma_t);
        REQUIRE(model.boundary_loss_correction == ground.boundary_loss_correction);
    }
    SECTION("empty cells inherit, given cells override, and the turbulence stays the run's")
    {
        const auto classes = ground_classes_from(read_ground_classes(folder.write(
                "class,name,impedance_model,flow_resistance,roughness\n"
                "0,grass,,,\n"
                "1,concrete,delany_bazley,10000,0\n")), ground, meters);
        REQUIRE(classes.size() == 2);
        REQUIRE(classes[0].model.flow_resistance == 200.0);
        REQUIRE(classes[0].model.impedance_model == ImpedanceModel::variable_porosity);
        REQUIRE(classes[0].model.roughness == 0.02);
        REQUIRE(classes[1].name == "concrete");
        REQUIRE(classes[1].line == 3);
        REQUIRE(classes[1].model.flow_resistance == 10000.0);
        REQUIRE(classes[1].model.impedance_model == ImpedanceModel::delany_bazley);
        REQUIRE(classes[1].model.roughness == 0.0);
        REQUIRE(classes[1].model.porosity_rate == 10.0);
        REQUIRE(classes[1].model.gamma_t == ground.gamma_t);
        REQUIRE(classes[1].model.receiver_height == ground.receiver_height);
        REQUIRE(find_ground_class(classes, 1) == &classes[1]);
        REQUIRE(find_ground_class(classes, 2) == nullptr);
    }
}

TEST_CASE("a ground class's lengths are in --units and nothing else is", "[ground_classes]")
{
    Folder folder;
    const GroundModel ground = run_ground();
    const auto file = read_ground_classes(folder.write("class,layer_depth,roughness,porosity_rate\n3,0.5,0.1,4\n"));
    const auto in_feet = ground_classes_from(file, ground, feet);
    // As the seam scales a command-line length: one multiply by the same constant.
    REQUIRE(in_feet[0].model.layer_depth == 0.5 * feet.to_meters);
    REQUIRE(in_feet[0].model.roughness == 0.1 * feet.to_meters);
    // A material rate, in 1/m whatever --units says.
    REQUIRE(in_feet[0].model.porosity_rate == 4.0);
    const auto in_meters = ground_classes_from(file, ground, meters);
    REQUIRE(in_meters[0].model.layer_depth == 0.5);
    REQUIRE(in_meters[0].model.roughness == 0.1);

    // A declared unit is checked against --units, and never converted.
    const auto declared = read_ground_classes(folder.write("class,roughness (ft)\n3,0.1\n"));
    REQUIRE(ground_classes_from(declared, ground, feet)[0].model.roughness == 0.1 * feet.to_meters);
    REQUIRE_THROWS_WITH(ground_classes_from(declared, ground, meters),
            ContainsSubstring("the roughness column declares ft, and this run reads lengths in m"));
    const auto depth = read_ground_classes(folder.write("class,layer_depth [m]\n3,0.1\n"));
    REQUIRE_THROWS_WITH(ground_classes_from(depth, ground, feet),
            ContainsSubstring("the layer_depth column declares m"));
}

TEST_CASE("a ground classes file refuses what it cannot read, naming the line", "[ground_classes]")
{
    Folder folder;
    const auto refuses = [&](const std::string& contents, const std::string& message) {
        INFO(contents);
        REQUIRE_THROWS_WITH(read_ground_classes(folder.write(contents)), ContainsSubstring(message));
    };
    refuses("class,flow\n0,200\n", "line 1: unknown column \"flow\". The columns are class, name");
    refuses("class,turbulence\n0,5e-4\n", "turbulence is the air's");
    refuses("class,Gamma_T\n0,1e-6\n", "turbulence is the air's");
    refuses("class,turbulence_model\n0,chessell\n", "turbulence is the air's");
    refuses("class,Class\n0,1\n", "column \"Class\" is given twice");
    refuses("name,flow_resistance\ngrass,200\n", "the header has no class column");
    refuses("class (ft),name\n0,grass\n", "column \"class (ft)\" declares ft, and it has no unit");
    refuses("class,flow_resistance (cgs)\n0,200\n", "declares cgs, and it is in kPa s/m^2");
    refuses("class,porosity_rate (1/ft)\n0,2\n", "declares 1/ft, and it is in 1/m");
    refuses("class,name\n\"0\",grass\n", "line 2 holds a quote");
    refuses("class,name\n0\n", "line 2 has 1 fields, and the header has 2");
    refuses("class,name\n0,grass,\n", "line 2 has 3 fields");
    refuses("class\n-1\n", "line 2: the class \"-1\" is not an integer from 0 to 2147483647");
    refuses("class\n1.0\n", "the class \"1.0\" is not an integer");
    refuses("class\n1e3\n", "the class \"1e3\" is not an integer");
    refuses("class\n\n#\n2147483648\n", "line 4: the class \"2147483648\" is not an integer");
    refuses("class\n\n", "has no classes");
    refuses("class,name\n0,grass\n,pad\n", "line 3: the class \"\" is not an integer");
    refuses("class,name\n0,grass\n1,grass\n0,pad\n", "the rows on lines 2 and 4 are both class 0");
    refuses("class,name\n0,my grass\n", "the name \"my grass\" may hold only");
    refuses("class,impedance_model\n0,delany-bazley\n", "line 2: unknown impedance model \"delany-bazley\"");
    refuses("class,flow_resistance\n0,0\n", "line 2: the flow resistance 0 kPa s/m^2 is not positive");
    refuses("class,flow_resistance\n0,-5\n", "is not positive");
    refuses("class,flow_resistance\n0,nan\n", "\"nan\", which is not a finite number");
    refuses("class,flow_resistance\n0,200 kPa\n", "\"200 kPa\", which is not a finite number");
    refuses("class,porosity_rate\n0,-1\n", "the porosity rate -1 is negative");
    refuses("class,layer_depth\n0,0\n", "the layer depth 0 is not positive");
    refuses("class,roughness\n0,-0.1\n", "the roughness -0.1 is negative");
    REQUIRE_THROWS_WITH(read_ground_classes(folder.write("")), ContainsSubstring("has no header"));
    REQUIRE_THROWS_WITH(read_ground_classes((folder.root / "absent.csv").string()),
            ContainsSubstring("Unable to open the ground classes file"));
}

TEST_CASE("a ground classes file keeps ids as land-cover codes, in file order", "[ground_classes]")
{
    Folder folder;
    // Ids are labels, not row indexes: any order, any gaps.
    const auto file = read_ground_classes(folder.write("class,name\n71,grassland\n21,developed\n2147483647,top\n"));
    REQUIRE(file.rows.size() == 3);
    REQUIRE(file.rows[0].id == 71);
    REQUIRE(file.rows[1].id == 21);
    REQUIRE(file.rows[2].id == 2147483647);
}

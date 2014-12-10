#include <vector>
#include <string>
#include <map>

typedef struct material
{
    int label;
    std::string name;
    double Er;
    double Ur;
    double econd;
    double mcond;

    material() {
        label = -1;
    };
    material(std::string name, double Er, double Ur, double econd, double mcond) : name(name), Er(Er), Ur(Ur), econd(econd), mcond(mcond) {
        label = -1;
    };
} MATERIAL;

typedef struct materiallist {
    std::map<std::string, int> labelmap;
    std::vector<MATERIAL> mats;

    void add_material(material m) {
        m.label = mats.size();
        labelmap[m.name] = mats.size();
        mats.push_back(m);
    }
} MATERIAL_LIST;

void add_material (MATERIAL_LIST &matl, const char *name, double Er, double Ur, double econd,double mcond);
void construct_mesh(MATERIAL_LIST &material_data, GRID* grid, SIMULATION* simulation);
int get_mlabel(MATERIAL_LIST &matl, const char *name);
void di_voxels(GRID* grid, SIMULATION* simulation, MATERIAL_LIST &mats, int is, int js, int ks, int iw, int jw, int kw, const char *);




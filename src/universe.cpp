#include "universe.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <omp.h>

#include "atom.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"
#include "reax_species.h"
#include "string_tools.h"

Universe::Universe(){};

Universe::~Universe() {
    prev_sys.reset();
    curr_sys.reset();
}

void Universe::flush() {
    prev_sys = curr_sys;
    curr_sys.reset();
}

void Universe::process_traj(const std::string &file_path, const std::vector<std::string> &type_names,
                            const float &rvdw_scale) {
    std::ifstream file(file_path);
    int curr_frame = 1;
    float neigh_radius = 2.5 * rvdw_scale;

    while (file.is_open() and !file.eof()) {
        flush();
        curr_sys = std::make_shared<System>();
        curr_sys->set_types(type_names);
        // curr_sys->link_moltopo_bank(moltopo_bank);

        if (file_path.ends_with(".lammpstrj"))
            curr_sys->load_lammpstrj(file);
        else if (file_path.ends_with(".xyz"))
            curr_sys->load_xyz(file);

        if (curr_sys->atoms.size() == 0) continue;

        // 使用 OpenMP 并行化计算密集型操作
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                curr_sys->search_neigh(neigh_radius, 10);
            }
            
            #pragma omp section
            {
                // 其他可以并行的预处理工作（如果有）
            }
        }
        
        // 这些操作有依赖关系，需要按顺序执行
        curr_sys->build_bonds_by_radius(rvdw_scale);
        curr_sys->build_molecules();

        // 分子公式计算可以并行
        std::vector<std::string> frame_formulas;
        frame_formulas.resize(curr_sys->molecules.size());
        
        #pragma omp parallel for
        for (size_t i = 0; i < curr_sys->molecules.size(); i++) {
            frame_formulas[i] = curr_sys->molecules[i]->formula;
        }
        
        reax_species->import_frame_formulas(frame_formulas);

        // 帧间依赖的操作，不能并行
        if (prev_sys && curr_sys) {
            update_reax_flow(prev_sys, curr_sys, curr_frame);
        }

        std::cout << "Frame: " << curr_frame << " ";
        curr_frame++;
        curr_sys->basic_info();
    }
    reax_species->analyze_frame_formulas();
}

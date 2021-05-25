#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;

int main (int argc, char* argv[]) {
   
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 100; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.4);

    auto start = std::chrono::steady_clock::now();
    auto X = bats::sample_sphere<double>(d, n);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Build data set needs: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;

	auto dist = bats::Euclidean();
    
    /*
    First we build a rips complex without extending filtrations
    */
	{
        start = std::chrono::steady_clock::now();
		auto R = bats::RipsComplex<CpxT>(X, dist, rmax, maxdim);
        std::cout << "\nRips complex summary:" << std::endl;
        end = std::chrono::steady_clock::now();
        auto t0 = end - start;
		R.print_summary();

        start = std::chrono::steady_clock::now();
		bats::ChainComplex<MT> C(R);
		auto RC = bats::ReducedChainComplex(C);        
        end = std::chrono::steady_clock::now();
        auto t1 = end - start;
        
        std::cout << "Reduced Chain Complex summary:" << std::endl;
        RC.print_summary();
        std::cout << "Build rips complex needs: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t0).count()
                << "ms" << std::endl;
        std::cout << "Build chain complex and reduce needs: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t1).count()
                << "ms" << std::endl;
	}

    /*
    Now let's exend a filtration, 
    for rips filtration, it is necessary to include all n choosing 2 edges to update filtration
    since we will not add/delete but only permute simplices for now.
    */
    start = std::chrono::steady_clock::now();
    //auto F = bats::RipsFiltration<CpxT>(X, dist, std::numeric_limits<double>::max(), maxdim);
    auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
    end = std::chrono::steady_clock::now();
    std::cout << "\nRips Filtration: with dimension" << std::endl;
    for (size_t i = 0; i <= F.maxdim(); i++) {
        std::cout << F.ncells(i) << " in dim " << i << std::endl;
    }
    std::cout << "Build rips filtration with infinite radius needs: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << "ms" << std::endl;
    
    /*
    Now let's reduce
    */
    start = std::chrono::steady_clock::now();
    auto FCC = bats::Chain(F, FT()); //FilteredChainComplex
    auto RFCC = bats::Reduce(FCC); //ReducedFilteredChainComplex    
    end = std::chrono::steady_clock::now();
    std::cout << "Build and Reduce Filtered Chain Complex"
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
            << "ms" << std::endl;

    std::cout << "Reduced matrix summary is" << std::endl;
    RFCC.print_summary();

    {
        /*
        If we want to change the filtration we have a variety of options.
        The naive way is to rebuild the filtration
        */
        // update filtration
        std::cout << "\nnaive way: rebuild the Filtered Chain Complex:" << std::endl;
        auto X_new = bats::sample_sphere<double>(d, n);
        auto F2 = bats::RipsFiltration<CpxT>(X_new, dist, std::numeric_limits<double>::max(), maxdim);
        
        start = std::chrono::steady_clock::now();
        auto FCC2 = bats::Chain(F2, FT()); //FilteredChainComplex
        auto RFCC2 = bats::Reduce(FCC2); //ReducedFilteredChainComplex
        end = std::chrono::steady_clock::now();
        std::cout << "Rebuild the filtration needs: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;

        std::cout << "Reduced matrix summary is" << std::endl;
        RFCC2.print_summary();
    }

    {   
        /*
        Second is to update the filtration on the FilteredChainComplex
        Then, we reduce the complex again
        */
        std::cout << "\nUpdate the filtration on the FilteredChainComplex" << std::endl;
        auto X_new = bats::sample_sphere<double>(d, n);
        auto F3 = bats::RipsFiltration<CpxT>(X_new, dist, std::numeric_limits<double>::max(), maxdim);
        auto vals = F3.vals();

        start = std::chrono::steady_clock::now();
        FCC.update_filtration(vals);
        auto RFCC3 = bats::Reduce(FCC);
        end = std::chrono::steady_clock::now();
        std::cout << "needs time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;
        
        std::cout << "Reduced matrix summary is" << std::endl;
        RFCC3.print_summary();

    }

    {
        /*
        The third option is to update the ReducedFilteredChainComplex
        In some situations, this may be faster
        */
        std::cout << "\nupdate the filtration on the ReducedFilteredChainComplex" << std::endl;
        auto X_new = bats::sample_sphere<double>(d, n);
        auto F4 = bats::RipsFiltration<CpxT>(X_new, dist, std::numeric_limits<double>::max(), maxdim);
        auto vals = F4.vals();

        start = std::chrono::steady_clock::now();
        RFCC.update_filtration(vals);
        end = std::chrono::steady_clock::now();
        std::cout << "needs time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;

        std::cout << "Reduced matrix summary is" << std::endl;
        RFCC.print_summary();
    }

    {
        /*
        The final option is to update the ReducedFilteredChainComplex using UQL
        In some situations, this may be faster
        */
        std::cout << "\nupdate the filtration on the ReducedFilteredChainComplex by UQL" << std::endl;
        auto X_new = bats::sample_sphere<double>(d, n);
        auto F5 = bats::RipsFiltration<CpxT>(X_new, dist, std::numeric_limits<double>::max(), maxdim);
        auto vals = F5.vals();

        start = std::chrono::steady_clock::now();
        RFCC.update_filtration_by_UQL(vals);
        end = std::chrono::steady_clock::now();
        std::cout << "needs time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
                << "ms" << std::endl;

        std::cout << "Reduced matrix summary is" << std::endl;
        RFCC.print_summary();
    }   

    return EXIT_SUCCESS;
}
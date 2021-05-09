#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>

#include <include/Matrix_type.hpp>

static constexpr int NI = 3;
static constexpr int NJ = 3;
static constexpr int NK = 3;
static constexpr int NL = 3;

static const std::array<size_t, 2> shape2D = { {NI, NJ} };
static const std::array<size_t, 3> shape3D = { {NI, NJ, NK} };
static const std::array<size_t, 4> shape4D = { {NI, NJ, NK, NL} };


int main()
{

    try {

        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_2D ..." << std::endl;
            MATRIX<2, int> A2(shape2D);
            std::cout << "CTOR MATRIX_2D ..." << std::endl;
            A2.Random(1, 8);
            std::cout << "RANDOM MATRIX_2D ..." << std::endl;
            std::cout << A2 << std::endl;
            std::cout << "DETERMINAT 2D = " << A2.DET() << std::endl;
            std::cout << "END TEST MATRIX_2D " << tmr.format() << std::endl;
        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_3D ..." << std::endl;
            MATRIX<3, int64_t> A3(shape3D);
            size_t count3D{};
            do {
                count3D++;
                A3.Random(1, 12);
                std::cout << "CTOR MATRIX_3D ..." << std::endl;
            } while (A3.DET_FULL() == 0);
            for(size_t i = 0; i < A3.size(0); ++i) {
                std::cout << "MATRIX_3D трансверсальные сечения ориентации i = \n" << A3.transversal_matrix('i', i) << std::endl;
            }
            for(size_t i = 1; i < A3.size(0); ++i) {
            std::cout << "Детерминант матрицы 3D ориентации i = " << A3.DET_orient('i', i) << std::endl;
            std::cout << "Детерминант матрицы 3D ориентации j = " << A3.DET_orient('j', i) << std::endl;
            std::cout << "Детерминант матрицы 3D ориентации k = " << A3.DET_orient('k', i) << std::endl;
            }
            std::cout << "Полный детерминант матрицы 3D = " << A3.DET_FULL() << std::endl;
            std::cout << "END TEST MATRIX_3D, ITERATION = " << count3D << ", " << tmr.format() << std::endl;

        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_4D ..." << std::endl;
            MATRIX<4, int64_t> A4(shape4D);
            size_t count4D{};
            do {
                count4D++;
                A4.Random(1, 12);
                std::cout << "CTOR MATRIX_4D ..." << std::endl;
            } while (A4.DET_FULL() == 0);
            for(size_t i = 0; i < A4.size(0); ++i) {
                std::cout << "MATRIX_4D трансверсальные сечения ориентации i = \n" << A4.transversal_matrix('i', i) << std::endl;
            }
            std::cout << "Вектор трансверсалей по альтернативному индексу i = ";
            for(size_t i = 1; i < A4.size(0); ++i)
                std::cout << A4.transversal_vector('i', i)[i] << " ";
            std::cout << std::endl;
            std::cout << "Вектор трансверсалей по альтернативному индексу j = ";
            for(size_t i = 1; i < A4.size(1); ++i)
                std::cout << A4.transversal_vector('j', i)[i] << " ";
            std::cout << std::endl;
            std::cout << "Вектор трансверсалей по альтернативному индексу k = ";
            for(size_t i = 1; i < A4.size(2); ++i)
                std::cout << A4.transversal_vector('k', i)[i] << " ";
            std::cout << std::endl;
            std::cout << "Вектор трансверсалей по альтернативному индексу l = ";
            for(size_t i = 1; i < A4.size(3); ++i)
                std::cout << A4.transversal_vector('l', i)[i] << " ";
            std::cout << std::endl;
            for(size_t i = 1; i < A4.size(0); ++i) {
                std::cout << "Детерминант матрицы 4D ориентации i = " << A4.DET_orient('i', i) << std::endl;
                std::cout << "Детерминант матрицы 4D ориентации j = " << A4.DET_orient('j', i) << std::endl;
                std::cout << "Детерминант матрицы 4D ориентации k = " << A4.DET_orient('k', i) << std::endl;
                std::cout << "Детерминант матрицы 4D ориентации l = " << A4.DET_orient('l', i) << std::endl;
            }

            std::cout << "Полный детерминант матрицы 4D = " << A4.DET_FULL() << std::endl;
            std::cout << "END TEST MATRIX_4D, ITERATION = " << count4D << ", " << tmr.format() << std::endl;
        }

    }
    catch (boost::system::system_error& ec) {
        std::cerr << "Error ocured! # " << ec.code() << " Messaga "
            << ec.what() << std::endl;
    }


    return 0;
}

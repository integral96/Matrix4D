#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>

//#include <include/GL_view.hpp>
#include <include/Matrix_type.hpp>

static constexpr int NI = 6;
static constexpr int NJ = 6;
static constexpr int NK = 6;
static constexpr int NL = 6;


static const std::array<size_t, 1> shape1D = {{NI}};
static const std::array<size_t, 2> shape2D = {{NI, NJ}};
static const std::array<size_t, 3> shape3D = {{NI, NJ, NK}};
static const std::array<size_t, 4> shape4D = {{NI, NJ, NK, NL}};

int main(int argc, char** argv)
{
//    {
//      glutInit(&argc, argv);
//      glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
//      glutInitWindowSize(windowWidth,windowHeight);
//      glutCreateWindow(windowName);
//      glClearColor(.9, .9, .9, 1);
//      glutDisplayFunc(display);
//      glutReshapeFunc(reshape);
//      glutKeyboardFunc(windowKey);
//      glutSpecialFunc(windowSpecial);

//      glutMainLoop();
//    }
    try {
        boost::timer::cpu_timer tmr;
        MATRIX<1, int> A1(shape1D);
        MATRIX<1, int> B1(shape1D);
        MATRIX<1, int> C1(shape1D);

        MATRIX<2, int> A2(shape2D);
        MATRIX<2, int> B2(shape2D);
        MATRIX<2, int> C2(shape2D);

        MATRIX<3, int> A3(shape3D);
        MATRIX<3, int> B3(shape3D);
        MATRIX<3, int> C3(shape3D);

        MATRIX<4, int> AA(shape4D);
        MATRIX<4, int> BB(shape4D);
        MATRIX<4, int> CC(shape4D);
        std::cout << "CTOR MATRIX" << tmr.format() << std::endl;
        A1.Random(0, 5);
        B1.Random(0, 3);

        A2.Random(0, 5);
        B2.Random(0, 3);

        A3.Random(0, 5);
        B3.Random(0, 3);

        AA.Random(0, 5);
        BB.Random(0, 3);
        std::cout << "RANDOM MATRIX" << tmr.format() << std::endl;
    //    std::cout << AA << std::endl;
    //    std::cout << BB << std::endl;
        std::cout << "START TEST MATRIX_1D ..." << std::endl;
        auto start1D = std::chrono::system_clock::now();
        size_t count1D{};
        do {
            B1 += (A1 + (B1*125));
            count1D++;
        } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start1D).count() < 20);
//        std::cout << B1 << std::endl;
        std::cout << "END TEST MATRIX_1D, ITERATION = " << count1D << tmr.format() << std::endl;
        std::cout << "START TEST MATRIX_2D ..." << std::endl;
        auto start2D = std::chrono::system_clock::now();
        size_t count2D{};
        do {
            C2 += ((A2*B2) + ((A2*B2)*125));
            count2D++;
        } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start2D).count() < 20);
//        std::cout << SD << std::endl;
        std::cout << "END TEST MATRIX_2D, ITERATION = " << count2D << tmr.format() << std::endl;
        std::cout << "START TEST MATRIX_3D ..." << std::endl;
        auto start3D = std::chrono::system_clock::now();
        size_t count3D{};
        do {
            C3 += (A3 + (B3*125));
            count3D++;
        } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start3D).count() < 20);
//        std::cout << SD << std::endl;
        std::cout << "END TEST MATRIX_3D, ITERATION = " << count3D << tmr.format() << std::endl;
        std::cout << "START TEST MATRIX_4D ..." << std::endl;
        auto start4D = std::chrono::system_clock::now();
        size_t count4D{};
        do {
            BB += (AA + (BB*125));
            count4D++;
        } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start4D).count() < 20);
//        std::cout << BB << std::endl;
        std::cout << "END TEST MATRIX_4D, ITERATION = " << count4D << tmr.format() << std::endl;
    //    for(const auto& x : BB.transversal_section('i').second) {
    //        std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << std::endl;
    //    }
    }  catch (boost::system::system_error& ec) {
        std::cerr << "Error ocured! # " << ec.code() << " Messaga "
                          << ec.what() << std::endl;
    }


    return 0;
}

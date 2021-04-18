#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

//#include <boost/timer/timer.hpp>

#include "GL_view.hpp"
//#include "Matrix4D.hpp"

//static constexpr int NI = 6;
//static constexpr int NJ = 6;
//static constexpr int NK = 6;
//static constexpr int NL = 6;

//using namespace _spatial;
//typedef boost::multi_array<int, 4> array_type;
//static const std::array<array_type::index, 4> shape = {{NI, NJ, NK, NL}};

int main(int argc, char** argv)
{
    {
      glutInit(&argc, argv);
      glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
      glutInitWindowSize(windowWidth,windowHeight);
      glutCreateWindow(windowName);
      glClearColor(.9, .9, .9, 1);
      glutDisplayFunc(display);
      glutReshapeFunc(reshape);
      glutKeyboardFunc(windowKey);
      glutSpecialFunc(windowSpecial);

      glutMainLoop();
    }
//    boost::timer::cpu_timer tmr;
//    Matrix4D<int> AA(shape);
//    Matrix4D<int> BB(shape);
//    Matrix4D<int> CC(shape);
//    std::cout << "BINGO_1" << tmr.format() << std::endl;
//    AA.Random(0, 5);
//    std::cout << "BINGO_1_ra" << tmr.format() << std::endl;
//    BB.Random(0, 3);
//    std::cout << "BINGO_2" << tmr.format() << std::endl;
////    std::cout << AA << std::endl;
////    std::cout << BB << std::endl;
//    auto start = std::chrono::system_clock::now();
//    do {
//        BB += (AA + (BB*125));
//    } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start).count() < 60);

//    std::cout << BB << std::endl;
////    for(const auto& x : BB.transversal_section('i').second) {
////        std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << std::endl;
////    }
//    std::cout << "BINGO_3" << tmr.format() << std::endl;

    return 0;
}

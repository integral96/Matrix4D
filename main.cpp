#include <iostream>
#include <vector>
#include <array>

#include "GL_view.hpp"

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

//    Matrix4D<int> AA(shape1);
//    Matrix4D<int> BB(shape1);
//    Matrix4D<int> CC(shape2);
//    AA.Random(0, 5);
//    BB.Random(0, 3);
//    BB += (AA - BB);
//    std::cout << BB << std::endl;
//    for(const auto& x : BB.transversal_section('i').second) {
//        std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << std::endl;
//    }

    return 0;
}

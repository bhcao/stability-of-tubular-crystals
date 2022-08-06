#include <cstdlib>
#include <fstream>

#include "figure.h"

void figure::draw(const char* name){
    std::ofstream fout("__nano__temp__.py");
    fout << "from matplotlib import pyplot as plt\n"
             "import numpy as np\nx = np.linspace(1, " 
         << this->data.size() << ", "
         << this->data.size() << ")\ny = np.array([";
    for (int i=0; i<this->data.size()-1; i++) {
        fout << this->data[i] << ", ";
    }
    fout << this->data[this->data.size()];
    fout << "])\nplt.figure()\n"
            "plt.plot(x, y)\nplt.savefig("
         << name << ")";
    fout.close();
    std::system("python __nano__temp__.py");
#ifdef _WIN32
    std::system("del __nano__temp__.py");
#endif
#ifdef __linux__
    std::system("rm __nano__temp__.py");
#endif
}
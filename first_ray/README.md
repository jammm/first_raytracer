# Compiling on OSX
(Warning: only works on jam's machine, You can fix it by changing the include and lib directories.)
`c++ -std=c++17 `pkg-config glfw3 assimp --cflags` *.cpp -o main `pkg-config glfw3 glew assimp --static --libs` -I./deps/cpp-taskflow-2.1.0  -framework OpenGL`

TODO: Use cmake

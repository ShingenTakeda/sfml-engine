#include "simple.hpp"

int main()
{
    Engine::Base app = {640, 640, "simple"};
    app.PushState(new SimpleState(&app));
    app.GameLoop();
}

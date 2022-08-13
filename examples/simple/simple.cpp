#include <iostream>
#include "simple.hpp"

void SimpleState::Init()
{
    std::cout << "Initializing state" << std::endl;
}

void SimpleState::Restart()
{
}


void SimpleState::CleanUp()
{
}

void SimpleState::Pause()
{
}

void SimpleState::Resume()
{
}

void SimpleState::Draw(double dt)
{
   std::cout << "Drawing" << std::endl;
}

void SimpleState::HandleInput(double dt)
{
     std::cout << "Handling input" << std::endl;
     if(this->base->b_event.type == sf::Event::Closed)
     {
         this->base->running = false;
     }
}

void SimpleState::Update(double dt)
{
    std::cout << "Updating" << std::endl;
}

SimpleState::SimpleState(Engine::Base *base)
{
    this->base = base;
}

SimpleState::~SimpleState()
{
    CleanUp();
}

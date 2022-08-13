#pragma once
#include <engine.hpp>

class SimpleState : public Engine::BState
{
    public:
        void Init() override;
        void CleanUp() override;
        void Restart() override;
        void Pause() override;
        void Resume() override;
        
        void Draw(double dt) override;
        void HandleInput(double dt) override;
        void Update(double dt) override;
        
        SimpleState(Engine::Base *base);
        ~SimpleState();
};

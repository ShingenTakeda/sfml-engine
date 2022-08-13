#pragma once
#include <iostream>
#include <lua.hpp>

namespace Engine
{
    //TODO:Look at elias dahler and LuaBridge for a good implementation of a Lua interface for OOP
    class ILua
    {
        public:
            ILua(const std::string &filename);
            ~ILua();
            
            void PrintError(const std::string &filename, const std::string& reason);
            
            template<typename T>
            T Get(const std::string &varName);
            
            bool GetToStack(const std::string &varName);
            
            template<typename T>
            T LuaGet(const std::string &varName) {return 0;};
            
            template<typename T>
            T LuaGetDefault(const std::string &varName) {return 0;};
            
            inline std::string LuaGetDefault() { return "null"; }
            
        private:
            lua_State *L;
            int level;
            std::string filename;
    };
};

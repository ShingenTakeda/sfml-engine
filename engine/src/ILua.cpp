#include <ILua.hpp>

Engine::ILua::ILua(const std::string &filename)
{
    L = luaL_newstate();
    if(luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0))
    {
        std::cout << "ERROR: no script loaded (" << filename << ")" << std::endl;
        L = nullptr;
    }
}

Engine::ILua::~ILua()
{
    if(L) lua_close(L);
}

void Engine::ILua::PrintError(const std::string &varName, const std::string &reason)
{
    std::cout << "ERROR: cannot get ["<< varName <<"]. " << reason << std::endl;
}

template<typename T> 
T Engine::ILua::Get(const std::string &varName)
{
    if(L == nullptr)
    {
        PrintError(varName, "No script loaded");
        return LuaGetDefault<T>(varName);
    }
    
    T result;
    if(GetToStack(varName))
    {
        result = LuaGet<T>(varName);
    }
    else
    {
        result = LuaGetDefault<T>();
    }
    
    lua_pop(L, level + 1);
    return result;
}

bool Engine::ILua::GetToStack(const std::string &varName)
{
    level = 0;
    std::string var = "";
        for(unsigned int i = 0; i < varName.size(); i++) {
      if(varName.at(i) == '.') {
        if(level == 0) {
          lua_getglobal(L, var.c_str());
        } else {
          lua_getfield(L, -1, var.c_str());
        }
 
        if(lua_isnil(L, -1)) {
          PrintError(varName, var + " is not defined");
          return false;
        } else {
          var = "";
          level++;
        }
      } else {
        var += varName.at(i);
      }
    }
    if(level == 0) {
      lua_getglobal(L, var.c_str());
    } else {
      lua_getfield(L, -1, var.c_str());
    }
    if(lua_isnil(L, -1)) {
        PrintError(varName, var + " is not defined");
        return false;
    }
 
    return true;
}

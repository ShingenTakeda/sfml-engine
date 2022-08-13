#pragma once
#include <map>
#include <vector>
#include <stack>
#include <cmath>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <lua.hpp>

namespace Engine
{
    /*
     * Base Engine stuff
     */
    class Base;
    
    const float gravityScale = 5.0f;
    const sf::Vector2f gravity( 0, 10.0f * gravityScale );
    const float dt = 1.0f / 60.0f;
    
    class BState
    {
        public:
            Base *base;
            
            virtual void Init() = 0;
            virtual void Restart(void) = 0;
            virtual void CleanUp() = 0;
            virtual void Pause() = 0;
            virtual void Resume() = 0;
            
            virtual void Draw(double dt) = 0;
            virtual void HandleInput(double dt) = 0;
            virtual void Update(double dt) = 0;
            
            BState() : base(nullptr) {}
            virtual ~BState() {}
    };
    
    class Base
    {
        public:
            std::stack<BState*> states;
            sf::RenderWindow window;
            sf::Event b_event;
            
            int screenWidth, screenHeight;
            
            bool running = true;
            bool vsync = true;
            double dt;
            
            void PushState(BState *state);
            void PopState();
            void ChangeState(BState *state);
            
            BState *PeekState();
            
            void GameLoop();
            
            Base(int w, int h, std::string wTitle);
            ~Base();
    };
    /*End*/
    
    /*
     * Dynamic input mapping
     */
    enum InputType
    {
      KEYBOARD,
      MOUSE,
      JOYSTICK
    };
    
    struct EKey
    {
        InputType inputType;
        sf::Event::EventType m_eventType;
        sf::Keyboard::Key m_keyCode;
        sf::Mouse::Button m_mouseButton;
    };
    
    extern std::map<std::string, EKey> map_keys;
    extern EKey key;
    /*End*/
    
    /*
     * Animation
     */
    enum AnimationMode
    {
      SINGLE,
      LOOP
    };
    
    struct Frame
    {
        sf::IntRect rect;
        double duration;
    };
    
    class Animation
    {
        private:
            std::vector<Frame> frames;
            double totalLength;
            double progress;
            sf::Sprite &target;
            
        public:
            Animation(sf::Sprite &target);
            virtual ~Animation() {}
            
            void AddFrame(Frame &frame);
            void Update(double elapsed);
            void Loop();
            const double GetLength() const {return totalLength;}
    };
    /*End*/
    
    /*
     * Math
     */
    
    extern const float PI;
    extern const float EPSILON;
    
#define EPSILON 0.0001f
#define PI 3.141592741f
    
#ifndef RAD2DEG
    float RAD2DEG(float radians);
#endif 
#ifndef DEG2RAD
    float DEG2RAD(float degrees);
#endif
    float CorrectDegrees(float degrees);
    
    inline float Random( float l, float h )
    {
        float a = (float)rand( );
        a /= (float)RAND_MAX;
        a = (h - l) * a + l;
        return a;
    }
    
    inline bool Equal(float a, float b)
    {
        return std::abs(a - b) <= EPSILON;
    }
    
    inline float Sqr(float a)
    {
        return a * a;
    }
    
    template <typename T>
    T DotProduct(const sf::Vector2<T> &l, const sf::Vector2<T> &r);
    
    template <typename T>
    T Magnitude(const sf::Vector2<T> &v);
    
    template <typename T>
    T Magnitude(const sf::Vector2<T> &a, const sf::Vector2<T> &b);
    
    template <typename T>
    T MagnitudeSqr(const sf::Vector2<T> &v);
    
    template <typename T>
    T MagnitudeSqr(const sf::Vector2<T> &a, const sf::Vector2<T> &b);
    
    template <typename T>
    void Normalize( sf::Vector2<T> &v);
    
    template <typename T>
    sf::Vector2<T> Normalized(const sf::Vector2<T> v);
    
    template <typename T>
    T Cross(const sf::Vector2<T> &l, const sf::Vector2<T> &r);
    
    template <typename T>
    sf::Vector2<T> Cross(const sf::Vector2<T> &v, T a);
    
    template <typename T>
    sf::Vector2<T> Cross(T a, const sf::Vector2<T> &v);
    
    template <typename T>
    T Angle(const sf::Vector2<T> &l, const sf::Vector2<T> &r);
    
    template <typename T>
    T Project(const sf::Vector2<T> &length, const sf::Vector2<T> &direction);
    
    template <typename T>
    T Perpendicular(const sf::Vector2<T> &len, const sf::Vector2<T> &dir);
    
    template <typename T>
    sf::Vector2<T> Reflection(const sf::Vector2<T> &vec, const sf::Vector2<T> &normal);
    
    /*End*/
    
    /*
     * Physics
     */
#define MAX_POLY_VERTEX_COUNT 64
    
    //TODO:Render stuff from a vertex array
    struct Body;
    
    template <typename T>
    struct TMat2
    {
        union
        {
            T v[4];
            T  m[2][2];
            
            struct
            {
                T right[2];
                T up[2];
            };
            
            struct
            {
                T xx; T xy;
                T yx; T yy;
            };
            
            struct
            {
                T c0r0; T c1r0;
                T c0r1; T c1r1;
            };
        };
        
        inline TMat2() :
        xx(1), xy(0),
        yx(0), yy(1) {}
        
        inline TMat2(T *fv) :
        xx(fv[0]), xy(fv[1]),
        yx(fv[2]), yy(fv[3]) {}
        
        inline TMat2
        (T _00, T _01,
         T _10, T _11) :
         xx(_00), xy(_01),
         yx(_10), yy(_11) {}
         
         void Set(float radians);
         
        const TMat2<T> operator*(const TMat2 &rhs ) const
        {
        // [00 01]  [00 01]
        // [10 11]  [10 11]

        return Mat2(
                v[0][0] * rhs[0][0] + v[0][1] * rhs[1][0],
                v[0][0] * rhs[0][1] + v[0][1] * rhs[1][1],
                v[1][0] * rhs[0][0] + v[1][1] * rhs[1][0],
                v[1][0] * rhs[0][1] + v[1][1] * rhs[1][1]);
        }
        
          const sf::Vector2<T> operator*( const sf::Vector2<T>& rhs ) const
        {
            return sf::Vector2<T>( xx * rhs.x + xy * rhs.y, yx * rhs.x + yy * rhs.y );
        }
    };

#define MAT2_EPSILON 0.000001f
    
    using Mat2 = TMat2<float>;
    using Mat2i = TMat2<int>;
    using Mat2d = TMat2<double>;
    
    struct Shape
    {
        enum Type
        {
            POLY,
            CIRCLE,
            COUNT
        };
        
        Shape() {}
        virtual Shape *Clone(void) const = 0;
        virtual void Initialize(void) = 0;
        virtual void ComputeMass(float density) = 0;
        virtual void SetOrientation(float rad) = 0;
        //virtual void DebugDraw(void) const = 0;
        virtual Type GetType(void) const = 0;
        
        Body *body;
        
        //For the circle shape
        float radius;
        
        //For polygon orientation
        Mat2 u;
    };
    
    struct Body
    {
        Body(Shape *shape_, uint32_t x, uint32_t y);
        
        void ApplyForce(const sf::Vector2f &f);
        void ApplyImpulse(const sf::Vector2f &impulse, const sf::Vector2f &contactVector);
        void SetStatic(void);
        void SetOrientation(float radians);
        
        sf::Vector2f position;
        sf::Vector2f velocity;
        
        float angularVelocity;
        float torque;
        float orientation;
        
        sf::Vector2f force;
        
        //Set by shape
        float I;//moment: inertia
        float iI;//inverse inertia
        float m;//mass
        float im;//inverse mass
        
        //For thes scene/jump table
        float staticFriction;
        float dynamicFriction;
        float restitution;
        
        Shape *shape;
        
        //For DebugDraw
        float r, g, b;
    };
    
    struct AABB
    {
        sf::Vector2f min;
        sf::Vector2f max;
    };
    
    struct Circle : public Shape
    {
        Circle(float r);
        
        Engine::Shape *Clone(void) const override;
        void Initialize(void) override;
        void ComputeMass(float density) override;
        void SetOrientation(float rad) override;
        //void DebugDraw(void) const override;
        Engine::Shape::Type GetType() const override;
        
    };
    
    struct PolygonShape : public Shape
    {
        void Initialize(void) override;
        Engine::Shape *Clone() const override;
        void ComputeMass(float density) override;
        void SetOrientation(float rad) override;
        Engine::Shape::Type GetType() const override;
        
        void SetBox(float hw, float hh);
        void Set(sf::Vector2f *vertices, uint32_t count);
        sf::Vector2f GetSupport(const sf::Vector2f &dir);
        
        uint32_t m_vertexCount;
        sf::Vector2f m_vertices[MAX_POLY_VERTEX_COUNT];
        sf::Vector2f m_normals[MAX_POLY_VERTEX_COUNT];
    };
    
    struct Manifold
    {
        Manifold(Body *_a, Body *_b) : a(_a), b(_b) {}
        
        void Solve(void);
        void Initialize(void);
        void ApplyImpulse(void);
        void PositionalCorrection(void);
        void InfiniteMassCorrection(void);
        
        Body *a;
        Body *b;
        
        float penetration;
        sf::Vector2f normal;
        sf::Vector2f contacts[2];
        uint32_t contact_count;
        float e;
        float df;
        float sf;
    };
    
    //Basic AABB stuff
    bool AABBvsAABB(AABB a, AABB b);
    
    //TODO; Implement  everything below here
    typedef void (*CollisionCallback)(Manifold *m, Body *a, Body *b);
    
    inline Engine::CollisionCallback Dispatch[Engine::Shape::COUNT][Engine::Shape::COUNT];
    
    //Maybe change this to make it accept only circles and polys
    void CircleCircle(Manifold *m, Body *a, Body *b);
    void CirclePolygon(Manifold *m, Body *a, Body *b);
    void PolygonCircle(Manifold *m, Body *a, Body *b);
    void PolygonPolygon(Manifold *m, Body *a, Body *b);
    
    //Physics data structure for physics simulation
    struct PWorld
    {
        PWorld(float _dt, uint32_t _iter) : p_dt(_dt), iterations(_iter) {}
        
        void Step(void);
        //void Render();
        Body *Add(Shape *shape, uint32_t x, uint32_t y);
        void Clear(void);
        
        float p_dt;
        uint32_t iterations;
        std::vector<Body *> bodies;
        std::vector<Manifold> contacts;
    };
    
    /*End*/
}

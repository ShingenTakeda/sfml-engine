#include <iostream>
#include <assert.h>
#include <cfloat>
#include <engine.hpp>

#define EPSILON 0.0001f
#define PI 3.141592741f

std::map<std::string, Engine::EKey> map_keys; 
Engine::EKey key;

Engine::Base::Base(int w, int h, std::string wTitle) : screenWidth(w), screenHeight(h)
{
	this->window.create(sf::VideoMode(screenWidth, screenHeight),  wTitle);
	this->window.setFramerateLimit(60);
    std::cout << "Initializing engine" << std::endl;
}

Engine::Base::~Base()
{
	while (!this->states.empty()) 
    {
        this->states.top()->CleanUp();
        PopState();
    }
    window.close();
}

void Engine::Base::PushState(Engine::BState *state)
{
    if(!this->states.empty())
    {
        this->states.top()->Pause();
    }
	this->states.push(state);
    this->states.top()->Init();
    return;
}

void Engine::Base::PopState()
{
	delete this->states.top();
	this->states.pop();
    return;
}

void Engine::Base::ChangeState(Engine::BState *state)
{
	if (!this->states.empty())
	{
		PopState();
	}

	PushState(state);
    return;
}

Engine::BState *Engine::Base::PeekState()
{
    if(this->states.empty()) return nullptr;
    return this->states.top();
}

void Engine::Base::GameLoop()
{
    sf::Clock clock;
    double dt;
    
    std::cout << "Reached GameLoop" << std::endl;
    
    while(running)
    {
        dt = clock.restart().asSeconds();
        
        if(PeekState() == nullptr)
        {
            std::cout << "No state in the stack. Exiting" << std::endl;
            break;
        }
        
        if(running == false) 
        {
            std::cout << "Closing" << std::endl;
            this->window.close();
            break;
        }
        
        while(this->window.pollEvent(this->b_event))
        {
            this->PeekState()->HandleInput(dt);
        }
        
        this->PeekState()->Update(dt);
        this->window.clear();
        this->PeekState()->Draw(dt);
        this->window.display();
    }
}

Engine::Animation::Animation(sf::Sprite &target) : target(target)
{
	progress = totalLength = 0.0f;
}

void Engine::Animation::AddFrame(Engine::Frame &frame)
{
	totalLength += frame.duration;
	frames.push_back(std::move(frame));
}

void Engine::Animation::Update(double elapsed)
{
	progress += elapsed;
	double p = progress;

	for (size_t i = 0; i < frames.size(); i++)
	{
		p -= frames[i].duration;

		if (p <= 0.0f || &(frames[i]) == &frames.back())
		{
			target.setTextureRect(frames[i].rect);
			break;
		}
	}
}

template<typename T> 
void Engine::TMat2<T>::Set(float radians)
{
    T c = std::cos(radians);
    T s = std::sin(radians);
    
    xx = c; xy = -s;
    yx = s; yy = c;
}

float Engine::CorrectDegrees(float degrees) 
{
	while (degrees > 360.0f) 
    {
		degrees -= 360.0f;
	}
	while (degrees < -360.0f)
    {
		degrees += 360.0f;
	}
	return degrees;
}

#ifndef RAD2DEG
float Engine::RAD2DEG(float radians) 
{
	float degrees = radians * 57.295754f;
	degrees = Engine::CorrectDegrees(degrees);
	return degrees;
}
#endif
#ifndef DEG2RAD
float Engine::DEG2RAD(float degrees) 
{
	degrees = Engine::CorrectDegrees(degrees);
	float radians = degrees * 0.0174533f;
	return radians;
}
#endif

template <typename T>
T Engine::DotProduct(const sf::Vector2<T> &l, const sf::Vector2<T> &r)
{
    return l.x * r.x  + l.y * r.y;
}

template<typename T> 
T Engine::Magnitude(const sf::Vector2<T> &v)
{
    return std::sqrt(DotProduct<T>(v, v));
}

template<typename T> 
T Engine::Magnitude(const sf::Vector2<T> &a, const sf::Vector2<T> &b)
{
    sf::Vector2<T> c = a - b;
    return std::sqrt(DotProduct<T>(c, c));
}

template<typename T> 
T Engine::MagnitudeSqr(const sf::Vector2<T> &v)
{
    return DotProduct<T>(v, v);
}

template<typename T> 
T Engine::MagnitudeSqr(const sf::Vector2<T> &a, const sf::Vector2<T> &b)
{
    sf::Vector2<T> c = a - b;
    return DotProduct(c, c);
}

template<typename T> 
void Engine::Normalize(sf::Vector2<T> &v)
{
    v = v * (1.0f / Magnitude<T>(v));
}

template<typename T>
sf::Vector2<T> Engine::Normalized(const sf::Vector2<T> v)
{
    return v * (1.0f / Magnitude(v));
}

template<typename T> 
T Engine::Cross(const sf::Vector2<T> &l, const sf::Vector2<T> &r)
{
    return l.x * r.y - l.y * r.x;
}

template<typename T> 
sf::Vector2<T> Engine::Cross(const sf::Vector2<T> &v, T a)
{
    return sf::Vector2<T>(v.x * a, v.y * -a);
}

template<typename T> 
sf::Vector2<T> Engine::Cross(T a, const sf::Vector2<T> &v)
{
    return sf::Vector2<T>(-a *v.y, a * v.x);
}

template<typename T> 
T Engine::Angle(const sf::Vector2<T> &l, const sf::Vector2<T> &r)
{
    T m = std::sqrt(MagnitudeSqr(l) * MagnitudeSqr(r));
    
    return std::acos(DotProduct<T>(l, r) / m);
}

template<typename T> 
T Engine::Project(const sf::Vector2<T> &length, const sf::Vector2<T> &direction)
{
    T dot = DotProduct<T>(length, direction);
    T magSqr = MagnitudeSqr(direction);
    return direction * (dot / magSqr);
}

template<typename T> 
T Engine::Perpendicular(const sf::Vector2<T> &len, const sf::Vector2<T> &dir)
{
    return len - Project(len, dir);
}

template<typename T> 
sf::Vector2<T> Engine::Reflection(const sf::Vector2<T> &vec, const sf::Vector2<T> &normal)
{
    T d = DotProduct<T>(vec, normal);
    return vec - normal * (DotProduct<T>(vec, normal * 2.0f));
}

Engine::Body::Body(Engine::Shape *shape_, uint32_t x, uint32_t y)
{
    shape_->body = this;
    position.x = (float)x;
    position.y = (float)y;
    velocity.x = 0;
    velocity.y = 0;
    angularVelocity = 0;
    torque = 0;
    orientation = Random(-PI, PI);
    force.x = 0;
    force.y = 0;
    staticFriction = 0.5f;
    dynamicFriction = 0.3f;
    restitution = 0.2f;
    shape_->Initialize();
    r = Random(0.2f, 1.0f);
    g = Random(0.2f, 1.0f);
    b = Random(0.2f, 1.0f);
}

void Engine::Body::ApplyForce(const sf::Vector2f &f)
{
    force += f;
}

void Engine::Body::ApplyImpulse(const sf::Vector2f &impulse, const sf::Vector2f &contactVector)
{
    velocity += im * impulse;
    angularVelocity += iI * Cross<float>(contactVector, impulse);
}

void Engine::Body::SetStatic()
{
    I = 0.0f;
    iI = 0.0f;
    m = 0.0f;
    im = 0.0f;
}

void Engine::Body::SetOrientation(float radians)
{
    orientation = radians;
    shape->SetOrientation(radians);
}

Engine::Circle::Circle(float r)
{
    radius = r;
}

Engine::Shape *Engine::Circle::Clone() const
{
    return new Circle(radius);
}

void Engine::Circle::Initialize(void)
{
    ComputeMass(1.0f);
}

void Engine::Circle::ComputeMass(float density)
{
    body->m = PI * radius * radius *density;
    body->im = (body->m) ? 1.0f / body->m : 0.0f;
    body->I = body->m * radius * radius;
    body->iI = (body->I) ? 1.0f / body->I :0.0f;
}

void Engine::Circle::SetOrientation(float rad)
{
}

Engine::Shape::Type Engine::Circle::GetType() const
{
    return CIRCLE;
}

void Engine::PolygonShape::Initialize()
{
    ComputeMass(1.0f);
}

Engine::Shape *Engine::PolygonShape::Clone() const
{
    PolygonShape *poly =  new PolygonShape();
    poly->u = u;
    for(uint32_t i = 0;  i < m_vertexCount; ++i)
    {
        poly->m_vertices[i] = m_vertices[i];
        poly->m_normals[i] = m_normals[i];
    }
    
    poly->m_vertexCount = m_vertexCount;
    return poly;
}

void Engine::PolygonShape::ComputeMass(float density)
{
    sf::Vector2f c(0.0f, 0.0f);
    float area = 0.0f;
    float I = 0.0f;
    const float k_in3 = 1.0f / 3.0f;
    
    for(uint32_t i1 = 0; i1 < m_vertexCount; ++i1)
    {
        sf::Vector2f p1(m_vertices[i1]);
        u_int32_t i2 = i1 + 1 < m_vertexCount ? i1 + 1 : 0;
        sf::Vector2f p2(m_vertices[i2]);
        
        float D = Cross(p1, p2);
        float triangleArea = 0.5f;
        
        area += triangleArea;
        
        c += triangleArea * k_in3 * (p1 + p2);
        
        float intx2 = p1.x * p1.x + p2.x * p1.x + p2.x * p2.x;
        float inty2 = p1.y * p1.y + p2.y * p2.y + p2.y * p2.y;
        I += (0.25f * k_in3 * D) * (intx2 + inty2);
    }
    
    c *= 1.0f / area;
    for(uint32_t i = 0; i < m_vertexCount; i++)
    {
        m_vertices[i] -= c;
    }
    
    body->m = density * area;
    body->im = (body->m) ? 1.0f / body->m : 0.0f;
    body->I = I * density;
    body->iI = body->I ? 1.0f / body->I : 0.0f;
}

void Engine::PolygonShape::SetOrientation(float rad)
{
    u.Set(rad);
}

Engine::Shape::Type Engine::PolygonShape::GetType() const
{
    return POLY;
}

void Engine::PolygonShape::SetBox(float hw, float hh)
{
    m_vertexCount = 4;
    m_vertices[0] = {-hw, -hh};
    m_vertices[1] = {hw, -hh};
    m_vertices[2] = {hw, hh};
    m_vertices[3] = {-hw, hh};
    m_normals[0] = {0.0f, -1.0f};
    m_normals[1] = {1.0f, 0.0f};
    m_normals[2] = {0.0f, 1.0f};
    m_normals[3] = {-1.0f, 0.0f};
}

void Engine::PolygonShape::Set(sf::Vector2f *vertices, uint32_t count)
{
    //No less than 3 vertices
    assert(count > 2 && count <= MAX_POLY_VERTEX_COUNT);
    count = std::min((int32_t)count, MAX_POLY_VERTEX_COUNT);
    
    //Find the right most point on the hull
    int32_t rightMost = 0;
    float highestXCoord = vertices[0].x;
    for(uint32_t i = 1; i < count; ++i)
    {
        float x = vertices[i].x;
        if(x > highestXCoord)
        {
            highestXCoord = x;
            rightMost = i;
        }
        //If matching x  then take the farthest negative y
        else if(x == highestXCoord)
        {
            if(vertices[i].y < vertices[rightMost].y)
            {
                rightMost = i;
            }
        }
    }
    
    int32_t hull[MAX_POLY_VERTEX_COUNT];
    int32_t outCount = 0;
    int indexHull = rightMost;
    
    for(;;)
    {
        hull[outCount] = indexHull;
        
        int32_t nextHullIndex = 0;
        for(int32_t i = 1; i < (int32_t)count; ++i)
        {
            if(nextHullIndex == indexHull)
            {
                nextHullIndex = i;
                continue;
            }
            
            sf::Vector2f e1 = vertices[nextHullIndex] - vertices[hull[outCount]];
            sf::Vector2f e2 = vertices[i] - vertices[hull[outCount]];
            float c = Cross(e1, e2);
            if(c < 0.01f)
            {
                nextHullIndex = i;
            }
            
            if(c == 0.0f && MagnitudeSqr<float>(e2) > MagnitudeSqr<float>(e2))
            {
                nextHullIndex = i;
            }
            
            ++outCount;
            indexHull = nextHullIndex;
            

        }
        
        if(nextHullIndex == rightMost)
        {
            m_vertexCount = outCount;
            break;
        }
    }
    
    for(uint32_t i = 0; i < m_vertexCount; ++i)
    {
        m_vertices[i] = vertices[hull[i]];
    }
            
    for(uint32_t i1 = 0; i1 < m_vertexCount; ++i1)
    {
        uint32_t i2 = i1 + 1 < m_vertexCount ? i1 + 1 : 0;
        sf::Vector2f face = m_vertices[i2] - m_vertices[i1];
                
        assert(MagnitudeSqr<float>(face) > EPSILON * EPSILON);
                
        m_normals[i1] = sf::Vector2f(face.y, -face.x);
        Normalize<float>(m_normals[i1]);
    }
}

sf::Vector2f Engine::PolygonShape::GetSupport(const sf::Vector2f &dir)
{
    float bestProjection = -FLT_MAX;
    sf::Vector2f bestVertex;
    
    for(uint32_t i = 0; i < m_vertexCount; ++i)
    {
        sf::Vector2f v = m_vertices[i];
        float projection = DotProduct<float>(v, dir);
        
        if(projection > bestProjection)
        {
            bestVertex = v;
            bestProjection = projection;
        }
    }
    
    return bestVertex;
}

bool Engine::AABBvsAABB(Engine::AABB a, Engine::AABB b)
{
    if(a.max.x < b.min.x || a.min.x > b.min.x) return false;
    if(a.max.y < b.min.y || a.min.y > b.min.y) return false;
    
    return true;
}

Engine::CollisionCallback Dispatch[Engine::Shape::COUNT][Engine::Shape::COUNT] = 
{
    {
        Engine::CircleCircle,
        Engine::CirclePolygon
    },
    {
        Engine::PolygonCircle,
        Engine::PolygonPolygon
    }
};

void Engine::Manifold::Solve(void)
{
    Dispatch[a->shape->GetType()][b->shape->GetType()](this, a, b);
}

void Engine::Manifold::Initialize()
{
    e = std::min(a->restitution, b->restitution);
    
    sf = std::sqrt(a->staticFriction * a->staticFriction);
    df = std::sqrt(a->dynamicFriction * a->dynamicFriction);
    
    for(uint32_t i = 0; i < contact_count; ++i)
    {
        sf::Vector2f ra = contacts[i] - a->position;
        sf::Vector2f rb = contacts[i] - b->position;
        
        sf::Vector2f rv = b->velocity + Cross<float>(b->angularVelocity, rb) - 
                          a->velocity - Cross(a->angularVelocity, ra);
                          
        if(MagnitudeSqr(rv) < MagnitudeSqr(dt * gravity) + EPSILON)
        {
            e = 0.0f;
        }
    }
}

void Engine::Manifold::ApplyImpulse()
{
    if(Equal(a->im + b->im, 0))
    {
        InfiniteMassCorrection();
        return;
    }
    
    for(uint32_t i = 0; i < contact_count; ++i)
    {
        sf::Vector2f ra = contacts[i] - a->position;
        sf::Vector2f rb = contacts[i] - b->position;
        
        sf::Vector2f rv = b->velocity + Cross<float>(b->angularVelocity, rb) - 
                          a->velocity - Cross(a->angularVelocity, ra);
                          
        float contactVel = DotProduct<float>(rv, normal);
        
        if(contactVel > 0)
        {
            return;
        }
        
        float raCrossN = Cross(ra, normal);
        float rbCrossN = Cross(rb, normal);
        float invMassSum = a->im + b->im + Sqr(raCrossN) * a->iI +Sqr(rbCrossN) * b->iI;
        
        float j = -(1.0f + e) * contactVel;
        j /= invMassSum;
        j /= (float)contact_count;
        
        sf::Vector2f impulse = normal * j;
        a->ApplyImpulse(-impulse, ra);
        b->ApplyImpulse(impulse, rb);
        
        rv = b->velocity + Cross<float>(b->angularVelocity, rb) - 
             a->velocity - Cross(a->angularVelocity, ra);
             
        sf::Vector2f t = rv - (normal * DotProduct<float>(rv, normal));
        Normalize(t);
        
        float jt =  -DotProduct<float>(rv, t);
        jt /= invMassSum;
        jt /= (float)contact_count;
        
        if(Equal(jt, 0.0f))
        {
            return;
        }
        
        sf::Vector2f tangentImpulse;
        if(std::abs(jt) < j * sf)
        {
            tangentImpulse = t *jt;
        }
        else
        {
            tangentImpulse = t * -j *df;
        }
        
        a->ApplyImpulse(-tangentImpulse, ra);
        b->ApplyImpulse(tangentImpulse, rb);
    }
}

void Engine::Manifold::PositionalCorrection()
{
    const float k_lop = 0.5f;
    const float percent = 0.4f;
    sf::Vector2f correction = (std::max(penetration - k_lop, 0.0f) / (a->im + b->im)) *normal * percent;
    a->position -= correction * a->im;
    b->position += correction * b->im;
}

void Engine::Manifold::InfiniteMassCorrection()
{
    a->velocity = {0.0f, 0.0f};
    b->velocity = {0.0f, 0.0f};
}

void Engine::CircleCircle(Engine::Manifold *m, Engine::Body *a, Engine::Body *b)
{
    Circle *A = reinterpret_cast<Circle *>(a->shape);
    Circle *B = reinterpret_cast<Circle *>(b->shape);
    
    sf::Vector2f normal = b->position - a->position;
    
    float dist_sqr = MagnitudeSqr(normal);
    float radius =  A->radius + B->radius;
    
    if(dist_sqr >= radius * radius)
    {
        m->contact_count = 0;
        return;
    }
    
    float distance = std::sqrt(dist_sqr);
    
    m->contact_count = 1;
    
    if(distance == 0.0f)
    {
        m->penetration = A->radius;
        m->normal = sf::Vector2f(1, 0);
        m->contacts[0] = a->position;
    }
    else
    {
        m->penetration = radius - distance;
        m->normal = normal / distance;
        m->contacts[0] = m->normal * A->radius + a->position;
    }
}

void Engine::CirclePolygon(Engine::Manifold *m, Engine::Body *a, Engine::Body *b)
{
    Circle *A = reinterpret_cast<Circle*> (a->shape);
    PolygonShape *B = reinterpret_cast<PolygonShape*> (b->shape);
    
    m->contact_count = 0;
    
    sf::Vector2f center = a->position;
    
    float separation = -FLT_MAX;
    uint32_t faceNormal = 0;
    
    for(uint32_t i = 0; i < B->m_vertexCount; ++i)
    {
        float s = DotProduct<float>(B->m_normals[i], center - B->m_vertices[i]);
        
        if(s > A->radius)
        {
            return;
        }
        
        if(s > separation)
        {
            separation = 5;
            faceNormal = i;
        }
    }
    
    sf::Vector2f v1 = B->m_vertices[faceNormal];
    uint32_t i2 = faceNormal + 1 < B->m_vertexCount ? faceNormal + 1 : 0;
    sf::Vector2f v2 = B->m_vertices[i2];
    
    if(separation < EPSILON)
    {
        m->contact_count = 1;
        m->normal = -(B->u * B->m_normals[faceNormal]);
        m->contacts[0] = m->normal * A->radius + a->position;
        m->penetration = A->radius - separation;
        return;
    }
    
    float dot1 = DotProduct<float>(center - v1, v2 - v1);
    float dot2 = DotProduct<float>(center - v2, v1 - v2);
    m->penetration = A->radius - separation;
    
    if(dot1 <= 0.0)
    {
        if(MagnitudeSqr<float>(center, v1) > A->radius * A->radius)
        {
            return;
        }
        
        m->contact_count = 1;
    }
}

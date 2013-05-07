// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif




////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

double RandomNumber(void)
{
#ifdef _WIN32
    // Seed random number generator
    static int first = 1;
    if (first) {
        srand(GetTickCount());
        first = 0;
    }
    
    // Return random number
    int r1 = rand();
    double r2 = ((double) rand()) / ((double) RAND_MAX);
    return (r1 + r2) / ((double) RAND_MAX);
#else
    // Seed random number generator
    static int first = 1;
    if (first) {
        struct timeval timevalue;
        gettimeofday(&timevalue, 0);
        srand48(timevalue.tv_usec);
        first = 0;
    }
    
    // Return random number
    return drand48();
#endif
}

////////////////////////////////////////////////////////////
// Finding Intersections
////////////////////////////////////////////////////////////
bool sphereIntersection(R3Sphere *sphere, R3Ray *ray, R3Intersection *intersection) {
    
    
    R3Vector L = sphere->Center() - ray->Start();
    double tca = L.Dot(ray->Vector());
    if (tca < 0)
        return false;
    
    double dSquared = (L.Length() * L.Length()) - (tca * tca);
    if (dSquared > (sphere->Radius() * sphere->Radius()))
        return false;
    
    double thc = sqrt((sphere->Radius() * sphere->Radius()) - dSquared);
    
    double t = min(tca + thc, tca - thc);
    
    if (t >= intersection->t)
        return false;
    
    
    intersection->t = t;
    intersection->hit = true;
    intersection->position = ray->Start() + (ray->Vector() * t);
    intersection->normal = intersection->position - sphere->Center();
    intersection->normal.Normalize();
    return true;
}

bool sphereIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    
    R3Sphere *sphere = node->shape->sphere;
    
    R3Vector L = sphere->Center() - ray->Start();
    double tca = L.Dot(ray->Vector());
    if (tca < 0)
        return false;
    
    double dSquared = (L.Length() * L.Length()) - (tca * tca);
    if (dSquared > (sphere->Radius() * sphere->Radius()))
        return false;
    
    double thc = sqrt((sphere->Radius() * sphere->Radius()) - dSquared);
    
    double t = min(tca + thc, tca - thc);
    
    if (t >= intersection->t)
        return false;
    
    
    intersection->t = t;
    intersection->hit = true;
    intersection->position = ray->Start() + (ray->Vector() * t);
    intersection->node = node;
    intersection->normal = intersection->position - sphere->Center();
    intersection->normal.Normalize();
    return true;
}

bool boxIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    R3Box *box = node->shape->box;
    bool updatedIntersection = false;
    
    //define planes
    vector<R3Plane *> planes;
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,0,1), box->Corner(0,1,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,1,0), box->Corner(1,0,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(1,0,0), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,1), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,0), box->Corner(0,1,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(1,1,0), box->Corner(1,0,0)));
    
    double epsilon = .000000001;
    
    for (int i = 0; i < (int)planes.size(); i++) {
        double denominator = ray->Vector().Dot(planes[i]->Normal());
        if (denominator == 0)
            continue;
        
        double numerator = planes[i]->Normal().Dot(planes[i]->Point() - ray->Start());
        double t = numerator/denominator;
        
        R3Point toCheck = ray->Start() + (ray->Vector() * t);
        
        if ((toCheck.X() <= box->XMin() - epsilon) || toCheck.X() > box->XMax() + epsilon)
            continue;
        
        if ((toCheck.Y() <= box->YMin() - epsilon) || toCheck.Y() > box->YMax() + epsilon)
            continue;
        
        if ((toCheck.Z() <= box->ZMin() - epsilon) || toCheck.Z() > box->ZMax() + epsilon)
            continue;
        
        if (t <= 0)
            continue;
        
        if (t < intersection->t) {
            intersection->t = t;
            intersection->hit = true;
            intersection->position = ray->Start() + (ray->Vector() * t);
            intersection->node = node;
            intersection->normal = planes[i]->Normal();
            intersection->normal.Normalize();
            updatedIntersection = true;
        }
    }
    
    for (int i = 0; i < (int)planes.size(); i++)
        delete planes[i];
    
    return updatedIntersection;
}

bool meshIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    R3Mesh *mesh = node->shape->mesh;
    bool updatedIntersection = false;
    double epsilon = .000001;
    
    for (int i = 0; i < (int)mesh->faces.size(); i++) {
        assert (mesh->faces[i]->vertices.size() == 3);
        
        
        //slight acceleration
        if (intersection->hit) {
            double dist = R3Distance(ray->Start(), mesh->faces[i]->vertices[0]->position);
            dist = min(dist, R3Distance(ray->Start(), mesh->faces[i]->vertices[1]->position));
            dist = min(dist, R3Distance(ray->Start(), mesh->faces[i]->vertices[2]->position));
            if (dist > intersection->t + epsilon)
                continue;
        }
        
        
        
        //check if ray intersects plane
        R3Plane *plane = &mesh->faces[i]->plane;
        double denominator = ray->Vector().Dot(plane->Normal());
        if (denominator == 0)
            continue;
        
        double numerator = plane->Normal().Dot(plane->Point() - ray->Start());
        double t = numerator/denominator;
        
        if (t <= 0)
            continue;
        
        //check if intersection is closer than previous found intersection
        if (intersection->t < t)
            continue;
        
        
        //check if intersects the face
        R3Vector v1 = mesh->faces[i]->vertices[0]->position - (ray->Start() + (ray->Vector() * t));
        R3Vector v2 = mesh->faces[i]->vertices[1]->position - (ray->Start() + (ray->Vector() * t));
        v1.Cross(v2);
        v1.Normalize();
        if (ray->Vector().Dot(v1) > 0)
            continue;
        
        v1 = mesh->faces[i]->vertices[1]->position - (ray->Start() + (ray->Vector() * t));
        v2 = mesh->faces[i]->vertices[2]->position - (ray->Start() + (ray->Vector() * t));
        v1.Cross(v2);
        v1.Normalize();
        if (ray->Vector().Dot(v1) > 0)
            continue;
        
        v1 = mesh->faces[i]->vertices[2]->position - (ray->Start() + (ray->Vector() * t));
        v2 = mesh->faces[i]->vertices[0]->position - (ray->Start() + (ray->Vector() * t));
        v1.Cross(v2);
        v1.Normalize();
        if (ray->Vector().Dot(v1) > 0)
            continue;
        
        
        //check if intersects the triangle
        
        
        intersection->t = t;
        intersection->face = mesh->faces[i];
        intersection->hit = true;
        intersection->position = ray->Start() + (ray->Vector() * t);
        intersection->node = node;
        intersection->normal = plane->Normal();
        intersection->normal.Normalize();
        updatedIntersection = true;
    }
    return updatedIntersection;
}

void faceIntersection(R3MeshFace *face, R3Ray *ray, R3Intersection *intersection, R3Node *lastNode) {
    
    if (!face)
        return;
    
    assert (face->vertices.size() == 3);
    
    
    //check if ray intersects plane
    R3Plane *plane = &face->plane;
    double denominator = ray->Vector().Dot(plane->Normal());
    if (denominator == 0)
        return;
    
    double numerator = plane->Normal().Dot(plane->Point() - ray->Start());
    double t = numerator/denominator;
    
    if (t <= 0)
        return;
    
    
    //check if intersection is closer than previous found intersection
    if (intersection->t < t)
        return;
    
    
    //check if intersects the face
    R3Vector v1 = face->vertices[0]->position - (ray->Start() + (ray->Vector() * t));
    R3Vector v2 = face->vertices[1]->position - (ray->Start() + (ray->Vector() * t));
    v1.Cross(v2);
    v1.Normalize();
    if (ray->Vector().Dot(v1) > 0)
        return;
    
    v1 = face->vertices[1]->position - (ray->Start() + (ray->Vector() * t));
    v2 = face->vertices[2]->position - (ray->Start() + (ray->Vector() * t));
    v1.Cross(v2);
    v1.Normalize();
    if (ray->Vector().Dot(v1) > 0)
        return;
    
    v1 = face->vertices[2]->position - (ray->Start() + (ray->Vector() * t));
    v2 = face->vertices[0]->position - (ray->Start() + (ray->Vector() * t));
    v1.Cross(v2);
    v1.Normalize();
    if (ray->Vector().Dot(v1) > 0)
        return;
    
    //check if intersects the triangle
    
    
    intersection->t = t;
    intersection->hit = true;
    intersection->position = ray->Start() + (ray->Vector() * t);
    intersection->node = lastNode;
    intersection->normal = plane->Normal();
    intersection->normal.Normalize();
}

bool cylinderIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    
    R3Cylinder *cylinder = node->shape->cylinder;
    bool updatedIntersection = false;
    
    //check mid section
    R3Point pA = cylinder->Axis().Start();
    R3Vector vA = cylinder->Axis().Vector();
    R3Point p = ray->Start();
    R3Vector v = ray->Vector();
    R3Vector deltaP = p - pA;
    
    
    double A = (v - v.Dot(vA) * vA).Dot(v - v.Dot(vA) * vA);
    double B = (v - v.Dot(vA) * vA).Dot(deltaP - deltaP.Dot(vA) * vA) * 2;
    double C = (deltaP - deltaP.Dot(vA) * vA).Dot(deltaP - deltaP.Dot(vA) * vA) - (cylinder->Radius() * cylinder->Radius());
    
    if ((B*B) - 4 * A * C >= 0) {
        double t1 = (-B + sqrt((B*B) - 4*A*C))/(2.*A);
        double t2 = (-B - sqrt((B*B) - 4*A*C))/(2.*A);
        
        if (t1 < 0)
            t1 = INFINITY;
        if (t2 < 0)
            t2 = INFINITY;
        double t = min(t1, t2);
        
        R3Point position = ray->Start() + (ray->Vector() * t);
        
        if (vA.Dot(position - cylinder->Axis().Start()) <= 0)
            t = INFINITY;
        
        if (vA.Dot(position - cylinder->Axis().End()) >= 0)
            t = INFINITY;
        
        
        //check top cap, start with plane intersection, see if within radius
        double t3 = INFINITY;
        bool intersectWithTopCap = false;
        double denominator = ray->Vector().Dot(cylinder->Axis().Vector());
        if (denominator != 0) {
            
            double numerator = cylinder->Axis().Vector().Dot(cylinder->Axis().Start() - ray->Start());
            double t3 = numerator/denominator;
            
            position = ray->Start() + (ray->Vector() * t3);
            
            if (R3Distance(position, cylinder->Axis().Start()) <= cylinder->Radius())
                if (t3 > 0) {
                    if (t3 < t) {
                        t = t3;
                        intersectWithTopCap = true;
                    }
                }
            
        }
        
        //check bottom cap, start with plane intersection, see if within radius
        bool intersectWithBottomCap = false;
        denominator = ray->Vector().Dot(cylinder->Axis().Vector());
        if (denominator != 0) {
            
            double numerator = cylinder->Axis().Vector().Dot(cylinder->Axis().End() - ray->Start());
            t3 = numerator/denominator;
            
            if (R3Distance(position, cylinder->Axis().End()) <= cylinder->Radius())
                if (t3 > 0) {
                    if (t3 < t) {
                        t = t3;
                        intersectWithBottomCap = true;
                    }
                }
        }
        
        
        if (t < intersection->t) {
            intersection->t = t;
            intersection->hit = true;
            intersection->position = ray->Start() + (ray->Vector() * t);
            intersection->node = node;
            updatedIntersection = true;
            if (intersectWithBottomCap)
                intersection->normal = cylinder->Axis().End() - cylinder->Center();
            else if (intersectWithTopCap)
                intersection->normal = cylinder->Axis().Start() - cylinder->Center();
            else {
                R3Point p1 = intersection->position;
                p1.Project(cylinder->Axis().Line());
                intersection->normal = intersection->position - p1;
            }
            intersection->normal.Normalize();
            
        }
    }
    return updatedIntersection;
}

bool boundingBoxIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    R3Box *box = &node->bbox;
    
    
    //define planes
    vector<R3Plane *> planes;
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,0,1), box->Corner(0,1,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,1,0), box->Corner(1,0,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(1,0,0), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,1), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,0), box->Corner(0,1,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(1,1,0), box->Corner(1,0,0)));
    
    double epsilon = .000000001;
    
    for (int i = 0; i < (int)planes.size(); i++) {
        double denominator = ray->Vector().Dot(planes[i]->Normal());
        if (denominator == 0)
            continue;
        
        double numerator = planes[i]->Normal().Dot(planes[i]->Point() - ray->Start());
        double t = numerator/denominator;
        
        if (t >= intersection->t)
            continue;
        
        R3Point toCheck = ray->Start() + (ray->Vector() * t);
        
        if ((toCheck.X() <= box->XMin() - epsilon) || toCheck.X() > box->XMax() + epsilon)
            continue;
        
        if ((toCheck.Y() <= box->YMin() - epsilon) || toCheck.Y() > box->YMax() + epsilon)
            continue;
        
        if ((toCheck.Z() <= box->ZMin() - epsilon) || toCheck.Z() > box->ZMax() + epsilon)
            continue;
        
        if (t <= 0)
            continue;
        
        if (t < intersection->t) {
            for (int x = 0; x < (int)planes.size(); x++)
                delete planes[x];
            
            return true;
        }
    }
    for (int x = 0; x < (int)planes.size(); x++)
        delete planes[x];
    return false;
}


bool ComputeIntersection(R3Node *node, R3Ray *ray, R3Intersection *intersection) {
    
    //recursive call
    if (node->children.size() > 0) {
        bool toReturn = false;
        for (int i = 0; i < (int)node->children.size(); i++) {
            ray->Transform(node->transformation.Inverse());
            bool recentReturn = false;
            if (boundingBoxIntersection(node->children[i], ray, intersection))
                recentReturn = ComputeIntersection(node->children[i], ray, intersection);
            ray->Transform(node->transformation);
            if (recentReturn) {
                intersection->position.Transform(node->transformation);
                intersection->normal.Transform(node->transformation.Inverse().Transpose());
                intersection->normal.Normalize();
                intersection->t = R3Distance(ray->Start(), intersection->position);
            }
            if (recentReturn)
                toReturn = recentReturn;
        }
        return toReturn;
    }
    
    
    bool updatedIntersection = false;
    if (node->shape != NULL) {
        switch(node->shape->type) {
            case R3_BOX_SHAPE:
                updatedIntersection = boxIntersection(node, ray, intersection);
                break;
            case R3_SPHERE_SHAPE:
                updatedIntersection = sphereIntersection(node, ray, intersection);
                break;
            case R3_CYLINDER_SHAPE:
                updatedIntersection = cylinderIntersection(node, ray, intersection);
                break;
            case R3_CONE_SHAPE:
                printf("Cone shape is not handled yet \n");
                break;
            case R3_MESH_SHAPE:
                updatedIntersection = meshIntersection(node, ray, intersection);
                break;
            case R3_SEGMENT_SHAPE:
                printf("Segment shape is not handled yet \n");
                break;
            default:
                return false;
                
        }
    }
    
    //transform intersection, normal, and parameter value t
    if (updatedIntersection) {
        intersection->position.Transform(node->transformation);
        intersection->normal.Transform(node->transformation.Inverse().Transpose());
        intersection->normal.Normalize();
        intersection->t = R3Distance(ray->Start(), intersection->position);
    }
    return updatedIntersection;
}



////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////
R3Particle* GenerateParticlesFromSphere(R3ParticleSource *pSource, double current_time, double delta_time) {
    
    assert (pSource->shape->type == R3_SPHERE_SHAPE);
    //generate random point on surface of sphere
    
    double u = (double)rand() / RAND_MAX;
    double v = (double)rand() / RAND_MAX;
    
    double theta = 2 * M_PI * u;
    double phi = acos((2 * v) - 1);
    double radius = pSource->shape->sphere->Radius();
    
    double x = pSource->shape->sphere->Center().X() + (radius * cos(theta) * sin(phi));
    double y = pSource->shape->sphere->Center().Y() + (radius * sin(theta) * sin(phi));
    double z = pSource->shape->sphere->Center().Z() + (radius * cos(phi));
    
    
    R3Particle *particle = new R3Particle();
    particle->position = R3Point(x, y, z);
    particle->velocity = pSource->velocity* (particle->position - pSource->shape->sphere->Center());
    particle->fixed = pSource->fixed;
    particle->drag = pSource->drag;
    particle->mass = pSource->mass;
    particle->material = pSource->material;
    particle->elasticity = pSource->elasticity;
    particle->lifetime = pSource->lifetime;
    pSource->numParticlesMade++; 
    return particle;
    
}

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
    
    // Generate new particles for every source
    for (int i = 0; i < scene->NParticleSources(); i++) {
        
        R3ParticleSource * pSource = scene->particle_sources[i];
        
        switch(pSource->shape->type) {
            case R3_BOX_SHAPE:
                break;
            case R3_SPHERE_SHAPE:
                while (pSource->numParticlesMade < pSource->rate * current_time)
                    scene->particles.push_back(GenerateParticlesFromSphere(pSource, current_time, delta_time));
                break;
            case R3_CYLINDER_SHAPE:
                break;
            case R3_CONE_SHAPE:
                break;
            case R3_MESH_SHAPE:
                break;
            case R3_SEGMENT_SHAPE:
                break;
            default:
                return;
                
        }
    }
}



////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////
void deleteParticle(R3Scene *scene, R3Particle *particle) {
    for (int i = 0; i < scene->NParticles(); i++) {
        if (scene->particles[i] == particle) {
            scene->particles[i] = scene->particles.back();
            scene->particles.pop_back();
            break;
        }
    }
}

R3Vector sphericalSink(R3Scene *scene, R3ParticleSink *sink, R3Particle *particle) {
    assert (sink->shape->type == R3_SPHERE_SHAPE);
    
    //fnd closest point
    R3Vector towardsSphere = sink->shape->sphere->Center() - particle->position;
    R3Ray *ray = new R3Ray(particle->position, towardsSphere);
    R3Intersection *intersection = new R3Intersection();
    intersection->t = INFINITY;
    sphereIntersection(sink->shape->sphere, ray, intersection);
    
    R3Vector towardsParticle = particle->position - sink->shape->sphere->Center();
    towardsParticle.Normalize();
    R3Point closestPoint = intersection->position;
    
    double d = R3Distance(closestPoint, particle->position);
    
    //delete particles that have collided with sink
    if (R3Distance(particle->position, sink->shape->sphere->Center()) < sink->shape->sphere->Radius())
        deleteParticle(scene, particle);
    
    
    //calculate force from sink
    towardsParticle.Flip();
    towardsParticle.Normalize();
    R3Vector force = towardsParticle;
    force *= sink->intensity/(sink->constant_attenuation + (sink->linear_attenuation * d) + (sink->quadratic_attenuation * d * d));
    force /= particle->mass;
    
    delete ray;
    delete intersection;
    
    return force;
}

R3Vector applyForces(R3Scene *scene, R3Particle *particle, double delta_time) {
    
    
    //apply gravity
    R3Vector force = scene->gravity * delta_time;
    
    //apply drag
    force -= (particle->velocity * particle->drag * delta_time)/(particle->mass);
    
    
    //apply force from sink
    for (int i = 0; i < scene->NParticleSinks(); i++) {
        
        R3ParticleSink *sink = scene->particle_sinks[i];
        
        switch(sink->shape->type) {
            case R3_BOX_SHAPE:
                break;
            case R3_SPHERE_SHAPE:
                force += sphericalSink(scene, sink, particle);
                break;
            case R3_CYLINDER_SHAPE:
                break;
            case R3_CONE_SHAPE:
                break;
            case R3_MESH_SHAPE:
                break;
            case R3_SEGMENT_SHAPE:
                break;
            default:
                return R3zero_vector;
                
        }
    }
    
    
    for (int i = 0; i < scene->NParticles(); i++) {
        if (particle == scene->particles[i])
            continue;
        R3Vector forceDirection = scene->particles[i]->position - particle->position;
        forceDirection.Normalize();
        double dist = R3Distance(scene->particles[i]->position, particle->position);
        double massProduct = scene->particles[i]->mass * particle->mass;
        double gravity = 6.67428e-11;
        if (dist == 0)
            continue;
        forceDirection *= (delta_time/particle->mass) * (gravity * massProduct)/(dist * dist);
        
        force += forceDirection;
        
    }
    
    //apply spring forces
    for (int i = 0; i < (int)particle->springs.size(); i++) {
        R3ParticleSpring *spring = particle->springs[i];
        R3Particle *p;
        R3Particle *q;
        
        //assign particles p and q
        if (spring->particles[0] == particle) {
            p = spring->particles[0];
            q = spring->particles[1];
        }
        else {
            p = spring->particles[1];
            q = spring->particles[0];
        }
        
   //     assert(!p->tail);
     //  assert(!q->tail);
  

        
        double d = R3Distance(p->position, q->position);
        if (d == 0)
            continue;
        R3Vector D = (q->position - p->position)/d;
        R3Vector velocDiff = q->velocity - p->velocity;
        
//            force += (delta_time/particle->mass) * ((spring->ks * (d - spring->rest_length))* D);
        force += (delta_time/particle->mass) * (((spring->ks * (d - spring->rest_length)) + (spring->kd * (velocDiff.Dot(D))))* D);
    }
    
    
    return force;
}

void UpdateParticle(R3Scene *scene, R3Particle *particle, double current_time, double delta_time, int integration_type) {
    
    double epsilon = .00000001;
    
    // Update position for every particle
    if (integration_type == EULER_INTEGRATION) {
        
        
        //check for intersection
        R3Ray *ray = new R3Ray(particle->position, particle->velocity);
        R3Intersection *intersection = new R3Intersection();
        intersection->t = INFINITY;
        ComputeIntersection(scene->root, ray, intersection);
        R3Point newPosition = particle->position + delta_time * particle->velocity;
        
        
        if (intersection->t < R3Distance(particle->position, newPosition)) {
            double elasped_time = R3Distance(intersection->position, particle->position)/particle->velocity.Length();
            assert (elasped_time < delta_time);
            //calculate specular velocity
            particle->position = intersection->position + (epsilon * intersection->normal);
            particle->velocity.Flip();
            particle->velocity.Rotate(intersection->normal, M_PI);
            //correct for elasticity
            R3Vector *temp = new R3Vector(particle->velocity);
            R3Vector verticalComp = *temp;
            R3Vector horizComp = *temp;
            verticalComp.Project(intersection->normal);
            horizComp -= verticalComp;
            verticalComp *= particle->elasticity;
            particle->velocity = horizComp + verticalComp;
            //calculate movement based on remaining time
            UpdateParticle(scene, particle, current_time + elasped_time, delta_time - elasped_time, integration_type);
            
        }
        
        // if no collision
        else {
            particle->position = newPosition;
            particle->velocity = particle->velocity + applyForces(scene, particle, delta_time);
            
        }
        delete ray;
        delete intersection;
    }
}

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{

    // Update position for every particle
    if (integration_type == EULER_INTEGRATION) {
        for (int i = 0; i < scene->NParticles(); i++) {

            //delete expired particles
            if ((scene->particles[i]->lifetime != 0) && (scene->particles[i]->lifetime <= delta_time))
                deleteParticle(scene, scene->particles[i]);
            else if (scene->particles[i]->lifetime != 0)
                scene->particles[i]->lifetime -= delta_time;
            

            if (!scene->particles[i]->fixed) {
                UpdateParticle(scene, scene->particles[i], current_time, delta_time, integration_type);

                
                R3Particle *tail = new R3Particle();
                tail->position = scene->particles[i]->position;
                tail->fixed = true;
                tail->tail = true;
                tail->lifetime = .5;
                tail->material = new R3Material();
                tail->material->kd[0] = scene->particles[i]->material->kd[0];
                tail->material->kd[1] = scene->particles[i]->material->kd[1];
                tail->material->kd[2] = scene->particles[i]->material->kd[2];
                tail->mass = 0;
                scene->particles.push_back(tail);
            }
        }
    }
}




////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time)
{
    // Draw every particle
    
    // REPLACE CODE HERE
    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int i = 0; i < scene->NParticles(); i++) {
        R3Particle *particle = scene->Particle(i);
        if (particle->tail)
            continue;
        glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
        const R3Point& position = particle->position;
        glVertex3d(position[0], position[1], position[2]);
    }
    glEnd();
}


void RenderTails(R3Scene *scene, double current_time, double delta_time)
{
    // Draw every particle
    
    // REPLACE CODE HERE
    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int i = 0; i < (int)scene->NParticles(); i++) {
        R3Particle *particle = scene->particles[i];
        if (!particle->tail)
            continue;
        double color0 = particle->material->kd[0] * particle->lifetime;
        double color1 = particle->material->kd[1] * particle->lifetime;
        double color2 = particle->material->kd[2] * particle->lifetime;
        glColor3d(color0, color1, color2);
        const R3Point& position = particle->position;
        glVertex3d(position[0], position[1], position[2]);
    }
    glEnd();
}




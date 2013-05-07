// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
using namespace std;
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif
#define PI 3.141592655359
#define GRAV_CONSTANT 6.67428e-11
#define ADAPTIVE_THRESHOLD 1e-2
#define eps 2e-12




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
// Generating Particles
////////////////////////////////////////////////////////////

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
    for (int i = 0; i < scene->NParticleSources(); i++) {
        
        R3ParticleSource *source = scene->ParticleSource(i);
        
        // spheres
        if (source->shape->type == R3_SPHERE_SHAPE) {
            
            // remainder term for fraction of particle
            double numberweshouldmake = delta_time * source->rate + source->remainder;
            int nparts = floor(numberweshouldmake);
            source->remainder = numberweshouldmake - nparts;
            
            for (int p = 0; p < nparts; p++) {
                R3Particle *newpart = new R3Particle();
                
                // calculate position and velocity
                double u = RandomNumber();
                double theta = 2*PI*u;
                double v = RandomNumber();
                double phi = acos(2*v-1);
                double r = source->shape->sphere->Radius();
                R3Point center = source->shape->sphere->Center();
                
                double x = r*cos(theta)*sin(phi);
                double y = r*sin(theta)*sin(phi);
                double z = r*cos(phi);
                
                R3Vector n = R3Vector(x,y,z);
                n.Normalize();
                
                // find tangent vector, use lecture notes to get velocity
                R3Plane plane = R3Plane(R3Point(0,0,0),n);
                R3Vector a;
                do {
                    a = R3Vector(RandomNumber(),RandomNumber(),RandomNumber());
                    a.Project(plane);
                } while (a.Length() == 0.0);
                a.Normalize();
                double t1 = 2*PI*RandomNumber();
                double t2 = sin(source->angle_cutoff)*RandomNumber();
                a.Rotate(n,t1);
                R3Vector vec = R3Vector(a);
                R3Vector cross = R3Vector(vec);
                cross.Cross(n);
                vec.Rotate(cross,acos(t2));
                
                newpart->position = center + n*r;
                newpart->velocity = vec*source->velocity;
                
                //update
                newpart->mass = source->mass;
                newpart->fixed = source->fixed;
                newpart->drag = source->drag;
                newpart->elasticity = source->elasticity;
                newpart->lifetime = source->lifetime;
                newpart->lifetimeactive = source->lifetimeactive;
                newpart->material = source->material;
                
                scene->particles.push_back(newpart);
            }
        }
        
        
        // CIRCLE
        if (source->shape->type == R3_CIRCLE_SHAPE) {
            
            double numberweshouldmake = delta_time * source->rate + source->remainder;
            int nparts = floor(numberweshouldmake);
            source->remainder = numberweshouldmake - nparts;
            
            for (int p = 0; p < nparts; p++) {
                R3Particle *newpart = new R3Particle();
                
                // calculate position and velocity
                
                double r = source->shape->circle->Radius();
                R3Point center = source->shape->circle->Center();
                R3Plane plane = source->shape->circle->Plane();
                R3Vector n = plane.Normal();
                n.Normalize();
                
                // get a random point on a circle
                double xcirc, ycirc;
                do {
                    xcirc = 2*r*(RandomNumber() - 0.5);
                    ycirc = 2*r*(RandomNumber() - 0.5);
                } while (xcirc*xcirc + ycirc*ycirc > r*r);
                
                // get basis vectors of circle
                R3Vector tang;
                do {
                    tang = R3Vector(RandomNumber(),RandomNumber(),RandomNumber());
                    tang.Project(plane);
                } while (tang.Length() == 0.0);
                R3Vector othertang = R3Vector(tang);
                othertang.Cross(n);
                othertang.Normalize();
                tang.Normalize();
                
                R3Point pos = center + tang*xcirc + othertang*ycirc;
                
                
                R3Vector a = R3Vector(tang);
                double t1 = 2*PI*RandomNumber();
                double t2 = sin(source->angle_cutoff)*RandomNumber();
                a.Rotate(n,t1);
                R3Vector vec = R3Vector(a);
                R3Vector cross = R3Vector(vec);
                cross.Cross(n);
                vec.Rotate(cross,acos(t2));
                
                newpart->position = pos;
                newpart->velocity = vec*source->velocity;
                
                
                newpart->mass = source->mass;
                newpart->fixed = source->fixed;
                newpart->drag = source->drag;
                newpart->elasticity = source->elasticity;
                newpart->lifetime = source->lifetime;
                newpart->lifetimeactive = source->lifetimeactive;
                newpart->material = source->material;
                
                scene->particles.push_back(newpart);
            }
        }
        
    }
    
}

////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

R3Vector ForceVector(R3Scene *scene, double current_time, R3Particle *particle, R3Point partpos, R3Vector partvel) {
    
    // set to 0
    R3Vector f = R3Vector(0,0,0);
    
    // gravity
    f += particle->mass*scene->gravity;
    
    // drag
    f += -1*partvel*particle->drag;
    
    // springs
    for (unsigned int i = 0; i < scene->particle_springs.size(); i++) {
        R3ParticleSpring *spring = scene->particle_springs[i];
        
        if (spring->particles[0] == particle && spring->particles[1] == particle) {
            
        }
        else if (spring->particles[0] == particle) {
            R3Point q = partpos;
            R3Point p = spring->particles[1]->position;
            R3Vector vq = partvel;
            R3Vector vp = spring->particles[1]->velocity;
            
            
            
            double d = R3Distance(p,q);
            if (d > eps) {
                R3Vector D = (q-p)/d;
                double mult = spring->kd*(vq - vp).Dot(D);
                f -= (spring->ks*(d - spring->rest_length) + mult)*D;
            }
        }
        
        else if (spring->particles[1] == particle) {
            R3Point q = partpos;
            R3Point p = spring->particles[0]->position;
            R3Vector vq = partvel;
            R3Vector vp = spring->particles[0]->velocity;
            double d = R3Distance(p,q);
            if (d > eps) {
                R3Vector D = (q-p)/d;
                double mult = spring->kd*(vq - vp).Dot(D);
                f -= (spring->ks*(d - spring->rest_length) + mult)*D;
            }
        }
    }
    
    // particle interactions
    
    for (int i = 0; i < scene->NParticles(); i++) {
        R3Particle *otherparticle = scene->Particle(i);
        if (otherparticle != particle) {
            R3Vector to = otherparticle->position - partpos;
            double d = to.Length();
            to.Normalize();
            if (d > eps) f += to*GRAV_CONSTANT*particle->mass*otherparticle->mass/d/d;
        }
    }
    
    // sink interaction
    for (int i = 0; i < scene->NParticleSinks(); i++) {
        R3ParticleSink *sink = scene->ParticleSink(i);
        
        //sphere sink
        if (sink->shape->type == R3_SPHERE_SHAPE) {
            // calculate closests point
            R3Point center = sink->shape->sphere->Center();
            double r = sink->shape->sphere->Radius();
            R3Vector to = partpos - center;
            to.Normalize();
            to *= r;
            R3Point closest = center + to;
            
            // calculate force
            R3Vector toward = closest - partpos;
            double d = toward.Length();
            toward.Normalize();
            if (d > eps) {
                toward *= sink->intensity/(sink->constant_attenuation + sink->linear_attenuation*d + sink->quadratic_attenuation*d*d);
                f += toward;
            }
        }
    }
    
    return f;
}



void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
    // array of ints to be deleted
    std::vector<int> tobedeleted;
    
    // for each particle
    for(int i = 0; i < scene->NParticles(); i++) {
        R3Particle *particle = scene->Particle(i);
        bool needstobedeleted = false;
        
        // check lifetime deletion
        if (particle->lifetimeactive) {
            if (particle->lifetime < 0) needstobedeleted = true;
            particle->lifetime -= delta_time;
        }
        
        R3Vector positionchange;
        R3Vector velocitychange;
        
        // Euler
        if (integration_type == EULER_INTEGRATION) {
            R3Vector f = ForceVector(scene, current_time, particle, particle->position, particle->velocity);
            positionchange = particle->velocity*delta_time;
            velocitychange = delta_time*f/particle->mass;
        }
        
        // midpoint
        else if (integration_type == MIDPOINT_INTEGRATION) {
            R3Vector f = ForceVector(scene, current_time, particle, particle->position, particle->velocity);
            R3Point xmid = particle->position + delta_time*particle->velocity/2.0;
            R3Vector vmid = particle->velocity + delta_time*f/particle->mass/2.0;
            R3Vector midf = ForceVector(scene, current_time, particle, xmid, vmid);
            positionchange = vmid*delta_time;
            velocitychange = delta_time*midf/particle->mass;
        }
        
        // rk4
        else if (integration_type == RK4_INTEGRATION) {
            R3Vector fk1 = ForceVector(scene, current_time, particle, particle->position, particle->velocity);
            R3Vector v1 = particle->velocity;
            R3Point x2 = particle->position + delta_time*particle->velocity/2.0;
            R3Vector v2 = particle->velocity + delta_time*fk1/particle->mass/2.0;
            R3Vector fk2 = ForceVector(scene, current_time, particle, x2, v2);
            R3Vector v3 = particle->velocity + delta_time*fk2/particle->mass/2.0;
            R3Vector fk3 = ForceVector(scene, current_time, particle, x2, v3);
            R3Point x4 = particle->position + delta_time*particle->velocity;
            R3Vector v4 = particle->velocity + delta_time*fk3/particle->mass;
            R3Vector fk4 = ForceVector(scene, current_time,particle,x4, v4);
            R3Vector finalf = (fk1+ 2*fk2 + 2*fk3 + fk4)/(6.0);
            R3Vector finalv = (v1+ 2*v2 +2*v3 + v4)/(6.0);
            
            positionchange = finalv*delta_time;
            velocitychange = delta_time*finalf/particle->mass;
        }
        
        
        // adaptive method
        else if (integration_type == ADAPTIVE_STEP_SIZE_INTEGRATION) {
            int numbersteps = 1;
            
            R3Vector poschange1;
            R3Vector velchange1;
            R3Vector poschange2;
            R3Vector velchange2;
            double error;
            
            // do first euler
            for (int step = 0; step < numbersteps; step++) {
                poschange2 = R3Vector(0,0,0);
                velchange2 = R3Vector(0,0,0);
                R3Vector f = ForceVector(scene, current_time, particle, particle->position + poschange2, particle->velocity + velchange2);
                poschange2 += particle->velocity*delta_time/numbersteps;
                velchange2 += delta_time*f/particle->mass/numbersteps;
            }
            // half below threshold
            do {
                poschange1 = poschange2;
                velchange1 = velchange2;
                numbersteps *= 2;
                
                // euler integrate with twice the number of steps
                for (int step = 0; step < numbersteps; step++) {
                    poschange2 = R3Vector(0,0,0);
                    velchange2 = R3Vector(0,0,0);
                    R3Vector f = ForceVector(scene, current_time, particle, particle->position + poschange2, particle->velocity + velchange2);
                    poschange2 += particle->velocity*delta_time/numbersteps/2.0;
                    velchange2 += delta_time*f/particle->mass/numbersteps/2.0;
                }
                
                double d1 = (poschange1 - poschange2).Length();
                double d2 = (velchange1 - velchange2).Length();
                
                // error for both velocity and position
                error = d1 + d2;
                
                
                
            } while (error > ADAPTIVE_THRESHOLD);
            
            positionchange = poschange2;
            velocitychange = velchange2;
            
            
        }
        
        else {
            fprintf(stderr, "invalid integration type\n");
            return;
        }
        
        
        // new values
        R3Point nextpos = particle->position + positionchange;
        R3Vector nextvel = particle->velocity + velocitychange;
        
        // check if particle needs to deleted
        for (int i = 0; i < scene->NParticleSinks(); i++) {
            R3ParticleSink *sink = scene->ParticleSink(i);
            
            //sphere sink
            if (sink->shape->type == R3_SPHERE_SHAPE) {
                double r = sink->shape->sphere->Radius();
                R3Point center = sink->shape->sphere->Center();
                
                // inside the sphere?
                if (R3Distance(nextpos,center) < r) needstobedeleted = true;
                
                R3Vector toward1 = particle->position - center;
                R3Vector toward2 = nextpos - center;
                
                // goes flying through the sphere?
                if (toward1.Dot(toward2) < 0) needstobedeleted = true;
                
            }
        }
        
        if (!needstobedeleted) {
            // check for bouncing in scene
            R3Vector to = nextpos - particle->position;
            to.Normalize();
            R3Ray *r = new R3Ray(particle->position, to);
            R3Intersect intersect = ComputeIntersect(scene,scene->Root(),r);
            
            if (intersect.intersected && intersect.t < R3Distance(nextpos,particle->position)) {
                
                R3Plane plane = R3Plane(intersect.pos,intersect.norm);
                
                R3Point first = particle->position;
                R3Point second = intersect.pos;
                R3Point third = nextpos;
                
                // calculate new position
                R3Vector posreflect = third-second;
                R3Vector posreflectplanar = R3Vector(posreflect);
                posreflectplanar.Project(plane);
                R3Vector posreflectnorm = posreflect - posreflectplanar;
                R3Vector purenorm = R3Vector(posreflectnorm);
                posreflectnorm *= -1*(particle->elasticity);
                purenorm.Normalize();
                purenorm *= -1;
                nextpos = second + posreflectplanar + posreflectnorm + eps*purenorm;
                
                // calculate new velocity
                
                R3Vector velreflect = R3Vector(nextvel);
                R3Vector velreflectplanar = R3Vector(velreflect);
                velreflectplanar.Project(plane);
                R3Vector velreflectnorm = velreflect - velreflectplanar;
                purenorm = R3Vector(velreflectnorm);
                velreflectnorm *= -1*(particle->elasticity);
                purenorm *= -1;
                purenorm.Normalize();
                nextvel = velreflectnorm + velreflectplanar + eps*purenorm;
                
            }
            
            
            // update position
            if (!particle->fixed) particle->position = nextpos;
            if (!particle->fixed) particle->velocity = nextvel;
            
            
        }
        
        // push on deletion vector
        else {
            int del = i;
            tobedeleted.push_back(del);
        }
        
    }
    
    // deleting particles
    for (unsigned int i = 0; i < tobedeleted.size(); i++) {
        //swap contents of partcile vector with last element
        R3Particle *temp = scene->particles.back();
        scene->particles[scene->particles.size() - 1] = scene->particles[tobedeleted[i]];
        scene->particles[tobedeleted[i]] = temp;
        
        // delete last element
        scene->particles.pop_back();
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
        glColor3d(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2]);
        const R3Point& position = particle->position;
        glVertex3d(position[0], position[1], position[2]);
    }   
    glEnd();
}




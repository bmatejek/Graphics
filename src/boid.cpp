//
//  boid.cpp
//
//
//  Created by Mark Whelan on 5/11/13.
//
//



#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
#include "boid.h"
#include <unistd.h>

#if defined(__APPLE__)
#define LINUX 0
#else
#define LINUX 1
#endif

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
// checking boid bullet intersections
////////////////////////////////////////////////////////////
double meshIntersection(R3Mesh *mesh, R3Ray *ray) {
    
    for (int i = 0; i < (int)mesh->faces.size(); i++) {
        assert (mesh->faces[i]->vertices.size() == 3);
        
        //check if ray intersects plane
        R3Plane *plane = &mesh->faces[i]->plane;
        double denominator = ray->Vector().Dot(plane->Normal());
        if (denominator == 0)
            continue;
        
        double numerator = plane->Normal().Dot(plane->Point() - ray->Start());
        double t = numerator/denominator;
        
        if (t <= 0)
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
        return t;
    }
    // return a ridiculously large number
    return 1e80;
}

double boxIntersection(R3Box bbox, R3Ray *ray) {
    R3Box *box = &bbox;

    double minIntersection = 1e80; 
    
    //define planes
    vector<R3Plane *> planes;
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,0,1), box->Corner(0,1,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(0,1,0), box->Corner(1,0,0)));
    planes.push_back(new R3Plane(box->Corner(0,0,0), box->Corner(1,0,0), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,1), box->Corner(0,0,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(0,1,0), box->Corner(0,1,1)));
    planes.push_back(new R3Plane(box->Corner(1,1,1), box->Corner(1,1,0), box->Corner(1,0,0)));
    
    double epsilon = 2e-12;
    
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
        
        if (t < minIntersection) {
            minIntersection = t; 
        }
    }
    
    for (int i = 0; i < (int)planes.size(); i++)
        delete planes[i];
    
    return minIntersection;
}

////////////////////////////////////////////////////////////
// Explode Boids
////////////////////////////////////////////////////////////
void deleteBoid(R3Scene *scene, R3Boid* boid) {
    
    for (int i = 0; i < (int)scene->boids.size(); i++) {
        if (scene->boids[i] == boid) {
            scene->boids[i] = scene->boids.back();
            scene->boids.pop_back();
            break;
        }
    }
}

void Explode(R3Scene *scene, R3Boid *boid) {
	if (boid->shape->type == R3_MESH_SHAPE) {
        int numParticlesperVertex = 2;
		for (unsigned int i = 0; i < boid->shape->mesh->vertices.size(); i++) {
            for (int x = 0; x < numParticlesperVertex; x++) {
                R3Particle *particle = new R3Particle();
                double speed = 1 * (double)rand() / RAND_MAX;;
                double x1 = 10 * (double)rand() / RAND_MAX;;
                double x2 = 10 * (double)rand() / RAND_MAX;;
                double x3 = 10 * (double)rand() / RAND_MAX;;
                double mass = 0.00000001;
                double drag = 0.0;
                double elasticity = 0.0;
                R3Vector velocity = R3Vector(x1, x2, x3);
                velocity.Normalize();
                
                static R3Material sink;
                static R3Material sink_material;
                static R3Material sink_material2;
                static R3Material sink_material3;
                
                if (sink.id != 33) {
                    sink.ka.Reset(0.2,0.2,0.2,1);
                    sink.kd.Reset(1,0,0,1);
                    sink.ks.Reset(1,0,0,1);
                    sink.kt.Reset(0,0,0,1);
                    sink.emission.Reset(1, 0, 0,1);
                    sink.shininess = 1;
                    sink.indexofrefraction = 1;
                    sink.texture = NULL;
                    sink.texture_index = -1;
                    sink.id = 33;
                }
                if (sink_material.id != 33) {
                    sink_material.ka.Reset(0.2,0.2,0.2,1);
                    sink_material.kd.Reset(1,0,0,1);
                    sink_material.ks.Reset(1,0,0,1);
                    sink_material.kt.Reset(0,0,0,1);
                    sink_material.emission.Reset(1, 0, 0,1);
                    sink_material.shininess = 1;
                    sink_material.indexofrefraction = 1;
                    sink_material.texture = NULL;
                    sink_material.texture_index = -1;
                    sink_material.id = 33;
                }
                if (sink_material2.id != 33) {
                    sink_material2.ka.Reset(0.2,0.2,0.2,1);
                    sink_material2.kd.Reset(0.96,0.44,0.11,1);
                    sink_material2.ks.Reset(0.96,0.44,0.11,1);
                    sink_material2.kt.Reset(0,0,0,1);
                    sink_material2.emission.Reset(0.96,0.44,0.11,1);
                    sink_material2.shininess = 1;
                    sink_material2.indexofrefraction = 1;
                    sink_material2.texture = NULL;
                    sink_material2.texture_index = -1;
                    sink_material2.id = 33;
                }
                if (sink_material3.id != 33) {
                    sink_material3.ka.Reset(0.2,0.2,0.2,1);
                    sink_material.kd.Reset(1,0.83,0,1);
                    sink_material.ks.Reset(1,0.83,0,1);
                    sink_material3.kt.Reset(0,0,0,1);
                    sink_material3.emission.Reset(1,0.83,0,1);
                    sink_material3.shininess = 1;
                    sink_material3.indexofrefraction = 1;
                    sink_material3.texture = NULL;
                    sink_material3.texture_index = -1;
                    sink_material3.id = 33;
                }
                
                particle->position = R3Point(boid->shape->mesh->vertices[i]->position);
                particle->velocity = speed * velocity;
                particle->mass = mass;
                particle->fixed = false;
                particle->drag = drag;
                particle->elasticity = elasticity;
                particle->lifetimeactive = true;
                particle->lifetime = 1.0;
                if (x1 < 0.5)
                    particle->material = &sink_material;
                else if (x1 < 3.33)
                    particle->material = &sink;
                else if (x1 < 6.67)
                    particle->material = &sink_material2;
                else
                    particle->material = &sink_material3;
                scene->particles.push_back(particle);
            }
	   }
	}
	pid_t pid;
	pid = fork();
	if (pid == 0) {
	  if (LINUX)
	    system("avplay -nodisp -autoexit boomBoid.wav");
	  else
	    system("afplay boomBoid.wav");
	  exit(0);
	}
}



////////////////////////////////////////////////////////////
// Generate Boids
////////////////////////////////////////////////////////////
void GenerateBoids(R3Scene *scene, int quantity, double distAway){
    
    
    for (int i = 0; i < quantity; i++) {
        double u = (double)rand() / RAND_MAX;
        double v = (double)rand() / RAND_MAX;
        
        double theta = 2 * M_PI * u;
        double phi = acos((2 * v) - 1);
        
        double x = scene->players[0]->pos.X() + (distAway * cos(theta) * sin(phi));
        double y = scene->players[0]->pos.Y() + (distAway * sin(theta) * sin(phi));
        double z = scene->players[0]->pos.Z() + (distAway * cos(phi));
  

        
        R3Boid *boid = new R3Boid();
        boid->pos = R3Point(x, y, z);
        
        //create mesh
        R3Mesh *mesh = new R3Mesh();
        if (!mesh) {
            fprintf(stderr, "Unable to allocate mesh\n");
            return;
        }
        
        // Read mesh file
        if (!mesh->Read("../input/shipAfirst.off")) {
            fprintf(stderr, "Unable to read mesh: ../input/shipAfirst.off\n");
            return;
        }
        
        
        // Create shape
        R3Shape *shape = new R3Shape();
        shape->type = R3_MESH_SHAPE;
        shape->box = NULL;
        shape->sphere = NULL;
        shape->cylinder = NULL;
        shape->cone = NULL;
        shape->mesh = mesh;
        shape->segment = NULL;
        boid->shape = shape;
        
        
        double error = .5;
        double v1 = (double)rand() / RAND_MAX;
        double mult1 = .6 + (v1 * error);
        
        //create boid
        R3Vector towardsPlayer = scene->players[0]->pos - boid->pos;
        towardsPlayer.Normalize();
        boid->velocity = towardsPlayer;
        boid->velocity.Normalize();
        boid->speed = mult1 * .6 * scene->players[0]->defaultVelocity;
        boid->health = 100;
        R3Material *material = new R3Material();
        material->kd[0] = 0;
        material->kd[1] = 0;
        material->kd[2] = 1;
        boid->material = material;
        
        //update mesh properties
        boid->shape->mesh->Translate(x, y, z);
        boid->shape->mesh->Center() = boid->pos;
        
        scene->boids.push_back(boid);
    }
}

//update boid 
void updateBoidVelocity(R3Scene *scene, R3Boid *boid) {
    
    double error = .6;
    double v1 = (double)rand() / RAND_MAX;
    double v2 = (double)rand() / RAND_MAX;
    double v3 = (double)rand() / RAND_MAX;
    double mult1 = 1 + (v1 * error);
    double mult2 = 1 + (v2 * error);
    double mult3 = 1 + (v3 * error);
    
    boid->velocity = scene->players[0]->pos - boid->pos;
    
    if (boid->velocity.Length() > 6) {
        
        boid->velocity.Normalize();
        
        boid->velocity.SetX(boid->velocity.X() * mult1);
        boid->velocity.SetY(boid->velocity.Y() * mult2);
        boid->velocity.SetZ(boid->velocity.Z() * mult3);
    }
    
    boid->velocity.Normalize();
}

//delete boids shot by player
void killShotBoids(R3Scene *scene, double delta_time) {
    double distAway = 15;
    
    for (int i = 0; i < (int)scene->bullets.size(); i++) {
        for (int j = 0; j < (int)scene->boids.size(); j++) {
            R3Ray *ray = new R3Ray(scene->bullets[i]->position, scene->bullets[i]->velocity);
            double intersection = meshIntersection(scene->boids[j]->shape->mesh, ray);
            if (intersection < scene->bullets[i]->velocity.Length() * delta_time) {
                Explode(scene, scene->boids[j]);
                deleteBoid(scene, scene->boids[j]);
                scene->players[0]->boidsKilled++;
                
                R3Bullet *temp = scene->bullets.back();
                scene->bullets[scene->bullets.size() - 1] = scene->bullets[i];
                scene->bullets[i] = temp;
                
                // delete last element
                scene->bullets.pop_back();
                if ((scene->players[0]->boidsKilled %5 == 0) && (scene->players[0]->boidsKilled != 0))
                    scene->players[0]->missiles++; 
                if (scene->boids.size() < 40)
                    GenerateBoids(scene, 2, distAway);
            }
        }
    }
}


bool ComputeBoidIntersection(R3Mesh *mesh, R3Player *player) {
    R3Box bbox = mesh->bbox;
    
    for (unsigned int i = 0; i < player->shape->mesh->vertices.size(); i++) {
        if (R3Distance(player->shape->mesh->vertices[i]->position, bbox.Centroid()) < bbox.DiagonalLength() / 2) {
            return true;
        }
    }
    
    return false;
}



void UpdateBoids(R3Scene *scene, double delta_time) {
    double distAway = 20;
    
    //check if any boids have been shot
    killShotBoids(scene, delta_time);
    
    //check if boid crashed into plane
    for (int i = 0; i < (int)scene->boids.size(); i++) {
        if (ComputeBoidIntersection(scene->boids[i]->shape->mesh, scene->players[0])) {
            scene->players[0]->health -= 5;
            Explode(scene, scene->boids[i]);
            deleteBoid(scene, scene->boids[i]);
            if (scene->boids.size() < 40)
                GenerateBoids(scene, 1, distAway);
        }
    }

    
    //update boid position and velocity
    for (int i = 0; i < (int)scene->boids.size(); i++) {
        scene->boids[i]->pos += delta_time * (scene->boids[i]->speed * scene->boids[i]->velocity);
        
        double dx = delta_time* scene->boids[i]->speed * scene->boids[i]->velocity.X();
        double dy = delta_time* scene->boids[i]->speed *scene->boids[i]->velocity.Y();
        double dz = delta_time* scene->boids[i]->speed * scene->boids[i]->velocity.Z();
        scene->boids[i]->shape->mesh->Translate(dx, dy, dz);
        scene->boids[i]->shape->mesh->Center() += delta_time * (scene->boids[i]->speed * scene->boids[i]->velocity);
        
        updateBoidVelocity(scene, scene->boids[i]);
    }
}

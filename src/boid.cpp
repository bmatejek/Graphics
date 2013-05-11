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

void GenerateBoids(R3Scene *scene, int quantity, double distAway){

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
    if (!mesh->Read("../input/smallTetra.off")) {
        fprintf(stderr, "Unable to read mesh: ../input/smallTetra.off\n");
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
    
    //create boid
    R3Vector towardsPlayer = scene->players[0]->pos - boid->pos;
    towardsPlayer.Normalize();
    boid->velocity = .5 * scene->players[0]->velocity * towardsPlayer;
    boid->health = 100;
    R3Material *material = new R3Material();
    material->kd[0] = 0;
    material->kd[1] = 0;
    material->kd[2] = 1;
    boid->material = material;
    
    scene->boids.push_back(boid);
    printf("num boids = %d \n", (int)scene->boids.size()); 
}

void UpdateBoids(R3Scene *scene, double delta_time) {
    
    for (int i = 0; i < (int)scene->boids.size(); i++) { 
        scene->boids[i]->pos += delta_time * (scene->boids[i]->velocity);
        
        double dx = delta_time* scene->boids[i]->velocity.X();
        double dy = delta_time* scene->boids[i]->velocity.Y();
        double dz = delta_time* scene->boids[i]->velocity.Z();
        scene->boids[i]->shape->mesh->Translate(dx, dy, dz);
        scene->boids[i]->shape->mesh->Center() += delta_time * (scene->boids[i]->velocity);
    }
}
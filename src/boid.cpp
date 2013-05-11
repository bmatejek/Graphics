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
    boid->position = R3Point(x, y, z);
    
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
    R3Vector towardsPlayer = scene->players[0]->pos - boid->position;
    towardsPlayer.Normalize();
    boid->velocity = .5 * scene->players[0]->velocity * towardsPlayer;
    boid->health = 100;
    R3Material *material = new R3Material();
    material->kd[0] = 0;
    material->kd[1] = 0;
    material->kd[2] = 1;
    boid->material = material;
    
    scene->boids.push_back(boid); 
}

void UpdateBoids(R3Scene *scene, double delta_time) {
    scene->players[0]->pos += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
    
    double dx = delta_time* scene->players[0]->velocity * scene->players[0]->nose.X();
    double dy = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Y();
    double dz = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Z();
    scene->players[0]->shape->mesh->Translate(dx, dy, dz);
    scene->players[0]->shape->mesh->Center() += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
}
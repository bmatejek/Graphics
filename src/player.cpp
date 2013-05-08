//
//  player.cpp
//  
//
//  Created by Mark Whelan on 5/7/13.
//
//

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
#include "player.h"
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

bool CheckBoundingBoxWithSphere(R3Sphere *sphere, R3Player *player) {
	R3Box bbox = player->shape->mesh->bbox;
	
	R3Point coords[8];
	
	coords[1] = bbox.Corner(0, 0, 0);
	coords[2] = bbox.Corner(0, 0, 1);
	coords[3] = bbox.Corner(0, 1, 0);
	coords[4] = bbox.Corner(0, 1, 1);
	coords[5] = bbox.Corner(1, 0, 0);
	coords[6] = bbox.Corner(1, 0, 1);
	coords[7] = bbox.Corner(1, 1, 0);
	coords[8] = bbox.Corner(1, 1, 1);
	for (int i = 0; i < 8; i++) {
		if (R3Distance(coords[i], sphere->Center()) < sphere->Radius())
			return true;
	}
	
	if (R3Distance(player->pos, sphere->Center()) < sphere->Radius())
		return true;
}

bool ComputeSphereIntersection(R3Sphere *sphere, R3Player *player) {
	if (CheckBoundingBoxWithSphere(sphere, player)) {
		for (unsigned int i = 0; i < player->shape->mesh->vertices.size(); i++) {
			if (R3Distance(player->shape->mesh->vertices[i]->position, sphere->Center()) < sphere->Radius())
				return true;
		}
	}
	return false;
}

bool ComputeMeshIntersection(R3Mesh *mesh, R3Player *player) {
	return false;
}

bool ComputeShapeIntersection(R3Scene *scene, R3Shape *shape, R3Player *player) {
	if (shape->type == R3_SPHERE_SHAPE) {
		return ComputeSphereIntersection(shape->sphere, player);
	}
	else if (shape->type == R3_MESH_SHAPE) {
		return ComputeMeshIntersection(shape->mesh, player);
	}
}

bool ComputeIntersection(R3Scene *scene, R3Node *node, R3Player *player) {
	// compute intersection for current object
	if (node->shape != NULL) {
		if (ComputeShapeIntersection(scene, node->shape, player))
			return true;
	}
	
	for (unsigned int i = 0; i < node->children.size(); i++) {
		if (ComputeIntersection(scene, node->children[i], player))
			return true;
	}
	
	return false;
}


void UpdatePlayers(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    scene->players[0]->pos += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
    
    double dx = delta_time* scene->players[0]->velocity * scene->players[0]->nose.X();
    double dy = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Y();
    double dz = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Z();
    scene->players[0]->shape->mesh->Translate(dx, dy, dz);
    scene->players[0]->shape->mesh->Center() += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
	
	if (ComputeIntersection(scene, scene->root, scene->players[0])) {
		printf("Here\n");
	}
}

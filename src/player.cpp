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
#include "bullet.h"
#include <sys/types.h>
#include <unistd.h>

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

double RandomNumber(void);


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

	// check if it is entirely within the sphere
	bool in = true;
	for (int i = 0; i < 8; i++) {
	  if (R3Distance(coords[i], sphere->Center()) < sphere->Radius()){
	    in = false;
	    break;
	  }
	}
	if (in) return false;

	for (int i = 0; i < 8; i++) {
		if (R3Distance(coords[i], sphere->Center()) < sphere->Radius())
			return true;
	}
	
	if (R3Distance(player->pos, sphere->Center()) < sphere->Radius())
		return true;
	return false;
}

bool ComputeSphereIntersection(R3Sphere *sphere, R3Player *player) {
	//if (CheckBoundingBoxWithSphere(sphere, player)) {
  
  bool in = true;
  bool out = true;
		for (unsigned int i = 0; i < player->shape->mesh->vertices.size(); i++) {
			
		  if (R3Distance(player->shape->mesh->vertices[i]->position, sphere->Center()) > sphere->Radius())
		    in = false;
		  if (R3Distance(player->shape->mesh->vertices[i]->position, sphere->Center()) < sphere->Radius())
		    out = false;
		}

		if (!in && !out) return true;
	//}
	return false;
}

bool ComputeMeshIntersection(R3Mesh *mesh, R3Player *player) {
	R3Box bbox = mesh->bbox;

	for (unsigned int i = 0; i < player->shape->mesh->vertices.size(); i++) {
		if (R3Distance(player->shape->mesh->vertices[i]->position, bbox.Centroid()) < bbox.DiagonalLength() / 2) {
			return true;
		}
	}
	
	return false;
}

bool ComputeShapeIntersection(R3Scene *scene, R3Shape *shape, R3Player *player) {
	if (shape->type == R3_SPHERE_SHAPE) {
		return ComputeSphereIntersection(shape->sphere, player);
	}
	else if (shape->type == R3_MESH_SHAPE) {
		return ComputeMeshIntersection(shape->mesh, player);
	}
	else {
		return false;
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

void Explode(R3Scene *scene, R3Player *player) {
	if (player->shape->type == R3_MESH_SHAPE) {
		for (unsigned int i = 0; i < player->shape->mesh->vertices.size(); i++) {
			int percent = 25;
			//int percent = player->shape->mesh->vertices.size() / 150;
			if (player->shape->mesh->vertices.size() < 300) {
				percent = 10;
			}
			if (i % percent == 0) {
				R3Particle *particle = new R3Particle();
				double speed = 1 * RandomNumber();
				double x1 = 10 * RandomNumber();
				double x2 = 10 * RandomNumber();
				double x3 = 10 * RandomNumber();
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
					sink_material2.emission.Reset(.96, .44, .11,1);
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
					sink_material3.emission.Reset(1, 0.83, 0,1);
					sink_material3.shininess = 1;
					sink_material3.indexofrefraction = 1;
					sink_material3.texture = NULL;
					sink_material3.texture_index = -1;
					sink_material3.id = 33;
				}
				
				particle->position = R3Point(player->shape->mesh->vertices[i]->position);
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
	  system("afplay boomPlayer.wav");
	  exit(0);
	}
	scene->players.erase(scene->players.begin());
}





void Explode(R3Scene *scene, R3Enemy *enemy) {
	if (enemy->shape->type == R3_MESH_SHAPE) {
		for (unsigned int i = 0; i < enemy->shape->mesh->vertices.size(); i++) {
			int percent = 1;
			//int percent = player->shape->mesh->vertices.size() / 150;
			/*if (enemy->shape->mesh->vertices.size() < 300) {
				percent = 10;
				}*/
			if (i % percent == 0) {
				R3Particle *particle = new R3Particle();
				double speed = 1 * RandomNumber();
				double x1 = 10 * RandomNumber();
				double x2 = 10 * RandomNumber();
				double x3 = 10 * RandomNumber();
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
					sink_material2.emission.Reset(.96, .44, .11,1);
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
					sink_material3.emission.Reset(1, 0.83, 0,1);
					sink_material3.shininess = 1;
					sink_material3.indexofrefraction = 1;
					sink_material3.texture = NULL;
					sink_material3.texture_index = -1;
					sink_material3.id = 33;
				}
				
				particle->position = R3Point(enemy->shape->mesh->vertices[i]->position);
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
	  system("afplay boomEnemy.wav");
	  exit(0);
	}

	scene->enemies.erase(scene->enemies.begin());
}







void UpdatePlayers(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    if (scene->players.size() != 0) {
		scene->players[0]->pos += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
		
		double dx = delta_time* scene->players[0]->velocity * scene->players[0]->nose.X();
		double dy = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Y();
		double dz = delta_time* scene->players[0]->velocity * scene->players[0]->nose.Z();
		scene->players[0]->shape->mesh->Translate(dx, dy, dz);
		scene->players[0]->shape->mesh->Center() += delta_time * (scene->players[0]->velocity * scene->players[0]->nose);
		
		if (ComputeIntersection(scene, scene->root, scene->players[0])) {
			Explode(scene, scene->players[0]);
		}
		
		for (unsigned int i = 0; i < scene->enemies.size(); i++) {
			// check collision with enemy fighters
			if (ComputeMeshIntersection(scene->enemies[i]->shape->mesh, scene->players[0])) {
				Explode(scene, scene->players[0]);
			}
		}
        
        if (scene->players[0]->missiletime > -1) {
            scene->players[0]->missiletime -= delta_time;
        }
        //fprintf(stdout,"%f %f %f %f %d\n",scene->players[0]->nose.X(),scene->players[0]->nose.Y(),scene->players[0]->nose.Z(),scene->players[0]->missiletime, scene->bullets.size());
    }
}

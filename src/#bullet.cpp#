//
//  bullet.cpp
//  
//
//  Created by Ethan Leeman on 5/7/13.
//
//

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "raytrace.h"
#include "bullet.h"
#include "particleview.h"
#include <cstdlib>
#include <time.h>
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
#define MISSILE_SCALE_FACTOR 0.1

static time_t last_bullet_sound = time(NULL);
static time_t last_missile_sound = time(NULL);
static bool bullet_shot = false;
static bool missile_shot = false;

void ShootBullet(R3Scene *scene) {
    //fprintf(stderr,"%d\n",scene->bullets.size());
    // generate a bullet from the plane
    
    
    R3Bullet *bullet;
    
    if (scene->players[0]->currentbullet == R3_REGULAR_BULLET) {
        bullet = new R3Bullet();
        bullet->type = scene->players[0]->currentbullet;
        bullet->position = scene->players[0]->pos + scene->players[0]->nose;
        bullet->velocity = 6*(scene->players[0]->velocity)*(scene->players[0]->nose);
        bullet->lifetimeactive = true;
        bullet->lifetime = 5.0;
        static R3Material sink_material;
        if (sink_material.id != 33) {
            sink_material.ka.Reset(0.2,0.2,0.2,1);
            sink_material.kd.Reset(1,0,0,1);
            sink_material.ks.Reset(1,0,0,1);
            sink_material.kt.Reset(0,0,0,1);
            sink_material.emission.Reset(0,0,0,1);
            sink_material.shininess = 1;
            sink_material.indexofrefraction = 1;
            sink_material.texture = NULL;
            sink_material.texture_index = -1;
            sink_material.id = 33;
        }
        bullet->material = &sink_material;

	// generate sound
	double seconds_since_start = difftime(time(NULL), last_bullet_sound);
	if (seconds_since_start > 0.5 || !bullet_shot) {
	  time(&last_bullet_sound);
	  bullet_shot = true;
	  pid_t pid;
	  pid = fork();
	  if (pid == 0) {
	    system("java BulletSound");
	    exit(0);
	  }
	}
    }
    
    if (scene->players[0]->currentbullet == R3_MISSILE_BULLET) {
        if (scene->players[0]->missiletime > 0) return;
        scene->players[0]->missiletime = 5.0;
        
        
        
        bullet = new R3Bullet();
        bullet->type = scene->players[0]->currentbullet;
        //create mesh
        R3Mesh *mesh = new R3Mesh();
        if (!mesh) {
            fprintf(stderr, "Unable to allocate mesh\n");
            return;
        }
        
        // Read mesh file
        if (!mesh->Read("../input/missile.off")) {
            fprintf(stderr, "Unable to read mesh: ../input/missile.off\n");
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
        
        bullet->shape = shape;
        bullet->position = scene->players[0]->pos + scene->players[0]->nose;
        bullet->velocity = 6*(scene->players[0]->velocity)*(scene->players[0]->nose);
        bullet->lifetimeactive = true;
        bullet->lifetime = 10.0;
        static R3Material sink_material;
        if (sink_material.id != 33) {
            sink_material.ka.Reset(0.2,0.2,0.2,1);
            sink_material.kd.Reset(1,0,0,1);
            sink_material.ks.Reset(1,0,0,1);
            sink_material.kt.Reset(0,0,0,1);
            sink_material.emission.Reset(0,0,0,1);
            sink_material.shininess = 1;
            sink_material.indexofrefraction = 1;
            sink_material.texture = NULL;
            sink_material.texture_index = -1;
            sink_material.id = 33;
        }
        bullet->material = &sink_material;
        bullet->shape->mesh->Scale(.3*MISSILE_SCALE_FACTOR,0.5*MISSILE_SCALE_FACTOR,MISSILE_SCALE_FACTOR);
        bullet->shape->mesh->Rotate(-PI/2,R3Line(R3Point(0,0,0), R3Point(0,0,1)));
        bullet->shape->mesh->Rotate(-PI/2,R3Line(R3Point(0,0,0), R3Point(1,0,0)));
        
        R3Vector begin = R3Vector(0,0,-1);
        R3Vector cross = scene->players[0]->nose;
        cross.Cross(begin);
        if (cross.Length() > 0) {
        R3Line line = R3Line(R3Point(0,0,0), cross);
        bullet->shape->mesh->Rotate(-acos(scene->players[0]->nose.Dot(begin)),line);
        }
        
        bullet->shape->mesh->Center() = bullet->position;
        double dx = bullet->position.X();
        double dy = bullet->position.Y();
        double dz = bullet->position.Z();
        bullet->shape->mesh->Translate(dx,dy,dz);
        
	double seconds_since_start = difftime(time(NULL), last_missile_sound);
	if (seconds_since_start > 4.0 || !missile_shot) {
	  time(&last_missile_sound);
	  missile_shot = true;
	  pid_t pid;
	  pid = fork();
	  if (pid == 0) {
	    system("java MissileSound");
	    exit(0);
	  }
	}
        
    }
    
    scene->bullets.push_back(bullet);
}

void UpdateBullets(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    std::vector<int> tobedeleted;
    int i;
    for (i = 0; i < (int)scene->bullets.size(); i++) {
        R3Bullet *bullet = scene->bullets[i];
    
        bullet->position += delta_time * bullet->velocity;
        if (bullet->type == R3_MISSILE_BULLET) {
            R3Vector change = delta_time * bullet->velocity;
            double dx = change.X();
            double dy = change.Y();
            double dz = change.Z();
            bullet->shape->mesh->Translate(dx,dy,dz);
        }
        
        if (bullet->lifetimeactive) {
            
            bullet->lifetime -= delta_time;
            
            if (bullet->lifetime < 0) {
                
                tobedeleted.push_back(i);
                
            }
            
        }
        
    }
    for (unsigned int i = 0; i < tobedeleted.size(); i++) {
        //swap contents of partcile vector with last element
        R3Bullet *temp = scene->bullets.back();
        scene->bullets[scene->bullets.size() - 1] = scene->bullets[tobedeleted[i]];
        scene->bullets[tobedeleted[i]] = temp;
        
        // delete last element
        scene->bullets.pop_back();
    }
}

/* Now in particleview
void RenderBullets(R3Scene *scene, double current_time, double delta_time)
{
    // Draw every particle
    
    // REPLACE CODE HERE
  //    glDisable(GL_LIGHTING);
    glPointSize(5);
    glBegin(GL_POINTS);

    
    for (int i = 0; i < (int)scene->bullets.size(); i += 10) {
        R3Bullet *bullet = scene->bullets[i];
	//    glColor3d(bullet->material->kd[0], bullet->material->kd[1], bullet->material->kd[2]);
	LoadMaterial(bullet->material);
        const R3Point& position = bullet->position;
        glVertex3d(position[0], position[1], position[2]);
    }
    glEnd();
}
*/

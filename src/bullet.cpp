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

void ShootBullet(R3Scene *scene) {
    //fprintf(stderr,"%d\n",scene->bullets.size());
    // generate a bullet from the plane
    R3Bullet *bullet = new R3Bullet();
    bullet->type = scene->players[0]->currentbullet;
    
    if (bullet->type == R3_REGULAR_BULLET) {
        bullet->position = scene->players[0]->pos + scene->players[0]->nose;
        bullet->velocity = 6*(scene->players[0]->velocity)*(scene->players[0]->nose);
        bullet->lifetimeactive = true;
        bullet->lifetime = 1.0;
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
    }
    
    if (bullet->type == R3_MISSILE_BULLET) {
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
        bullet->lifetime = 1.0;
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
    }
    
    scene->bullets.push_back(bullet);
}

void UpdateBullets(R3Scene *scene, double current_time, double delta_time, int integration_type) {
    
    std::vector<int> tobedeleted;
    int i;
    for (i = 0; i < (int)scene->bullets.size(); i++) {
        R3Bullet *bullet = scene->bullets[i];
    
        bullet->position += delta_time * bullet->velocity;
        
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

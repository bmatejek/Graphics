// Source file for raytracing code


// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"
#include <limits> 
using namespace std;



////////////////////////////////////////////////////////////////////////
// Create image from scene
//
// This is the main ray tracing function called from raypro
// 
// "width" and "height" indicate the size of the ray traced image
//   (keep these small during debugging to speed up your code development cycle)
//
// "max_depth" indicates the maximum number of secondary reflections/transmissions to trace for any ray
//   (i.e., stop tracing a ray if it has already been reflected max_depth times -- 
//   0 means direct illumination, 1 means one bounce, etc.)
//
// "num_primary_rays_per_pixel" indicates the number of random rays to generate within 
//   each pixel during antialiasing.  This argument can be ignored if antialiasing is not implemented.
//
// "num_distributed_rays_per_intersection" indicates the number of secondary rays to generate
//   for each surface intersection if distributed ray tracing is implemented.  
//   It can be ignored otherwise.
// 
////////////////////////////////////////////////////////////////////////

//quadratic formula, first answer, -1 if no answer
static double Quad1(double a, double b, double c) {
    double d = b*b-4*a*c;
    if (d < 0) return -1;
    double total = (-b + sqrt(d))/(2*a);
    return total;
    
}

//quadratic formula, second answer, -1 if no answer
static double Quad2(double a, double b, double c) {
    double d = b*b-4*a*c;
    if (d < 0) return -1;
    double total = (-b - sqrt(d))/(2*a);
    return total;
}


// create rays through a pixel, as in lecture
R3Ray ConstructRayThroughPixel(int x, int y, int width, int height, R3Point eye, R3Vector towards, R3Vector up, R3Vector right, double xfov, double yfov) {
    
    R3Point p = towards.Point();
    
    p.Rotate(up,xfov*(1-(2.0*x)/width));
    p.Rotate(right,-1*(yfov*(1-(2.0*y)/height)));

    return R3Ray(eye,eye+p);
}

// Create the specular ray, by mirroring over tangent plane
R3Ray SpecularRay(R3Ray r, R3Intersect hit) {
    R3Point p = hit.pos;
    R3Vector vect = R3Vector(r.Vector());
    R3Plane plane = R3Plane(hit.norm,0);
    vect.Mirror(plane);
    return R3Ray(p, vect);
}

// Create Transmission ray, by simply changing the hit position
R3Ray TransmissionRay(R3Ray r, R3Intersect hit) {
    R3Point p = hit.pos;
    R3Vector vect = R3Vector(r.Vector());
    return R3Ray(p, vect);
}

// Create Refraction Ray using Snell's law
R3Ray RefractionRay(R3Ray r, R3Intersect hit, double indexinc, double indexout) {
    R3Vector n = hit.norm;
    R3Vector l = r.Vector();
    n.Normalize();
    l.Normalize();
    l *= -1;
    double cosi = n.Dot(l);
    double sini = sqrt(1-cosi*cosi);
    double sinr = indexinc*sini / indexout;
    double cosr = sqrt(1-sinr*sinr);
    R3Vector t = (indexinc*cosi/indexout - cosr)*n - indexinc*l/indexout;
    
    R3Point p = hit.pos;
    t.Normalize();
    return R3Ray(p, t);
}

// find intersection with plane, as in lecture
R3Intersect PlaneIntersect(R3Ray r, R3Plane p) {
    R3Intersect sect = R3Intersect();
    
    R3Point p0 = r.Start();
    R3Vector V = r.Vector();
    R3Vector N = p.Normal();
    double d = p.D();
    
    // parallel lines
    if (V.Dot(N) == 0) {
        sect.intersected = false;
        return sect;
    }
    double t = -1*(N.Dot(p0.Vector())+d)/(V.Dot(N));
    R3Vector posv = p0.Vector()+t*V;
    R3Point pos = posv.Point();
    
    sect.intersected = true;
    sect.pos = pos;
    sect.norm = N;
    
    if (sect.norm.Dot(r.Vector()) > 0) sect.norm *= -1;
    
    sect.t = t;
    
    return sect;
}


// intersection with cylinder
R3Intersect CylinderIntersect(R3Ray r, R3Cylinder cyl) {
    
    // circular part
    R3Intersect sect = R3Intersect();
    sect.intersected = false;
    sect.t = -1;
    
    double x0 = r.Start().X();
    double z0 = r.Start().Z();
    double xc = cyl.Center().X();
    double yc = cyl.Center().Y();
    double zc = cyl.Center().Z();
    double xr = r.Vector().X();
    double zr = r.Vector().Z();
    double rad = cyl.Radius();
    
    double a = zr*zr+xr*xr;
    double b = 2*xr*(x0-xc)+2*zr*(z0-zc);
    double c = (x0-xc)*(x0-xc)+(z0-zc)*(z0-zc)-rad*rad;
    double t1 = Quad1(a,b,c);
    double t2 = Quad2(a,b,c);
    
    double t = min(t1,t2);

    // if one answer is negative, either ray starts inside, or never hits
    if (t < 0) {
        if (t1 < 0 && t2 < 0) return sect;
        
        t = max(t1,t2);

    }
    R3Point pt = r.Point(t);
    
    if (yc - rad/2 < pt.Y() && pt.Y() < yc + rad/2) {
        sect.intersected = true;
        sect.t = t;
        sect.pos = pt;
        R3Vector v = pt - cyl.Center();
        v.SetY(0);
        v.Normalize();
        sect.norm = v;
    }
    
    R3Intersect secttop = R3Intersect();
    secttop.intersected = false;
    secttop.t = -1;
    
    
    // top face
    R3Point ptop = cyl.Center() + R3Point(0,cyl.Height()/2,0);
    R3Vector normal = R3Vector(0,1,0);
    R3Plane pltop = R3Plane(ptop,normal);
    secttop = PlaneIntersect(r,pltop);
    double d = R3Distance(secttop.pos, ptop);
    if (d > cyl.Radius()) {
        secttop.intersected = false;
    }
    
    R3Intersect sectbottom = R3Intersect();
    sectbottom.intersected = false;
    sectbottom.t = -1;
    
    
    // bottom face
    R3Point pbottom = cyl.Center() + R3Point(0,-cyl.Height()/2,0);
    normal = R3Vector(0,-1,0);
    R3Plane plbottom = R3Plane(pbottom,normal);
    sectbottom = PlaneIntersect(r,plbottom);
    d = R3Distance(sectbottom.pos, pbottom);
    if (d > cyl.Radius()) {
        sectbottom.intersected = false;
    }
    
    if (secttop.intersected && !sect.intersected) sect = secttop;
    if (secttop.intersected && sect.intersected) {
        if (secttop.t < sect.t) sect = secttop;
    }
    if (sectbottom.intersected && !sect.intersected) sect = sectbottom;
    if (sectbottom.intersected && sect.intersected) {
        if (sectbottom.t < sect.t) sect = sectbottom;
    }
    
    return sect;
}

// intersect with cone
R3Intersect ConeIntersect(R3Ray r, R3Cone cone) {
    R3Intersect sect = R3Intersect();
    sect.intersected = false;
    sect.t = -1;
    double eps = 1e-12;
    
    double h = cone.Height();
    double px = cone.Center().X();
    double py = cone.Center().Y() + h/2;
    double pz = cone.Center().Z();
    double rx = r.Start().X();
    double ry = r.Start().Y();
    double rz = r.Start().Z();
    double vx = r.Vector().X();
    double vy = r.Vector().Y();
    double vz = r.Vector().Z();
    double rad = cone.Radius();
    
    
    double a = vx*vx/rad/rad + vz*vz/rad/rad-vy*vy/h/h;
    double b = 2*(px-rx)*(-1*vx)/rad/rad + 2*(pz-rz)*(-1*vz)/rad/rad - 2*(py-ry)*(-1*vy)/h/h;
    double c = (px-rx)*(px-rx)/rad/rad +(pz-rz)*(pz-rz)/rad/rad -(py-ry)*(py-ry)/h/h;
    double t1 = Quad1(a,b,c);
    double t2 = Quad2(a,b,c);
    
    double t = min(t1,t2);
    
    if (t < -eps) {
     if (t1 < -eps && t2 < -eps) return sect;
        t = max(t1,t2);
    }

    R3Point pt = r.Point(t);
    
    if (pt.Y() > py || pt.Y() < py - h) {
        sect.intersected = false;
    }
    
    else {
    sect.intersected = true;
    sect.t = t;
    sect.pos = pt;
    R3Vector v = pt - cone.Center();
    v.SetY(0);
    v.Normalize();
    sect.norm = v;
    }
    
    
    R3Intersect basesect = R3Intersect();
    basesect.intersected = false;
    basesect.t = -1;
    
    R3Point midbase = R3Point(px,py-h,pz);
    R3Plane plane = R3Plane(midbase,R3Vector(0,-1,0));
    basesect = PlaneIntersect(r,plane);
    if (R3Distance(basesect.pos,midbase) > rad) {
        basesect.intersected = false;
    }
    
    if (basesect.intersected && sect.intersected) {
        if (basesect.t < sect.t) return basesect;
        return sect;
    }
    if (basesect.intersected) {
        return basesect;
    }
    return sect;
    
    
    return sect;
}

// intersect with box
R3Intersect BoxIntersect(R3Ray r, R3Box b) {
    double eps = 1e-12;
    
    R3Plane *p = new R3Plane[6];
    p[0] = R3Plane(b.Min(), R3Vector(0,0,-1));
    p[1] = R3Plane(b.Min(), R3Vector(0,-1,0));
    p[2] = R3Plane(b.Min(), R3Vector(-1,0,0));
    p[3] = R3Plane(b.Max(), R3Vector(0,0,1));
    p[4] = R3Plane(b.Max(), R3Vector(0,1,0));
    p[5] = R3Plane(b.Max(), R3Vector(1,0,0));
    R3Intersect answer = R3Intersect();
    answer.intersected = false;
    double tmin = std::numeric_limits<double>::infinity();
    R3Intersect *i = new R3Intersect[6];
    for (int j = 0; j < 6; j++) {
        i[j] = PlaneIntersect(r,p[j]);
        
        if (i[j].intersected == true) {
            if (i[j].t < tmin && i[j].t > 0) {
                
                if (j % 3 == 0) {
                    if (i[j].pos.X() > b.Min().X() -eps && i[j].pos.X() < b.Max().X() +eps ) {
                        if (i[j].pos.Y() > b.Min().Y() -eps && i[j].pos.Y() < b.Max().Y() +eps)
                        {
                            answer = i[j];
                            tmin = i[j].t;
                        }
                    }
                }
                else if (j % 3 == 1) {
                    if (i[j].pos.X() > b.Min().X() -eps && i[j].pos.X() < b.Max().X()+eps ) {
                        if (i[j].pos.Z() > b.Min().Z()-eps && i[j].pos.Z() < b.Max().Z()+eps )
                        {
                            answer = i[j];
                            tmin = i[j].t;
                        }
                    }
                }
                else if (j % 3 == 2) {
                    if (i[j].pos.Z() > b.Min().Z()-eps && i[j].pos.Z() < b.Max().Z()+eps ) {
                        if (i[j].pos.Y() > b.Min().Y()-eps && i[j].pos.Y() < b.Max().Y() +eps)
                        {
                            answer = i[j];
                            tmin = i[j].t;
                        }
                    }
                }
            }
        }
    }
    return answer;
}

// intersect with a triangle, geometrically as in lecture
R3Intersect TriIntersect(R3Ray r, R3MeshFace f) {
    
    
    R3Plane plane = f.plane;
    R3Intersect sect = PlaneIntersect(r,plane);
    if (sect.intersected == false) return sect;
    
    for (int i = 0; i < 3; i++) {
        
        R3Vector p0 = r.Start().Vector();
        R3Vector v1 = f.vertices[i]->position.Vector() + (-1)*p0;
        R3Vector v2 = f.vertices[(i+1)%3]->position.Vector() + (-1)*p0;
        v2.Cross(v1);
        R3Vector n1 = v2;
        n1.Normalize();
        
        
        R3Plane p = R3Plane(r.Start(),n1);
        if (R3SignedDistance(p,sect.pos) < 0) { sect.intersected = false; return sect; }
        
    }
    return sect;
    
}

// intersect with a triangle mesh, by going through all the triangles
R3Intersect TriMeshIntersect(R3Ray r, R3Mesh m) {
    double currentmin = std::numeric_limits<double>::infinity();
    R3Intersect finalsect;
    finalsect.intersected = false;
    finalsect.t = currentmin;
    for (int i = 0; i < m.NFaces(); i++) {
        R3MeshFace *f = m.Face(i);
        R3Intersect sect = TriIntersect(r,*f);
       
        if (sect.intersected == true) {
            
            if (sect.t < currentmin && sect.t > 0) {
                currentmin = sect.t;
                finalsect = sect;
                
            }
        }
    }
    return finalsect;
}

// intersect with a sphere, as in lecture
R3Intersect SphereInt(R3Ray r, R3Sphere s) {
    double eps = 1e-12;
    R3Intersect i = R3Intersect();
    double d = R3Distance(s.Center(),r);

    if (d > s.Radius()) {
        i.intersected = false;
    }
    
    else if (R3Distance(r.Start()+eps*r.Vector(),s.Center()) < s.Radius()) {
        i.intersected = true;
        
        
        double xc = r.Start().X() - s.Center().X();
        double yc = r.Start().Y() - s.Center().Y();
        double zc = r.Start().Z() - s.Center().Z();
        double xt = r.Vector().X();
        double yt = r.Vector().Y();
        double zt = r.Vector().Z();
        
        double a = xt*xt + yt*yt + zt*zt;
        double b = 2*xt*xc + 2*yt*yc + 2*zt*zc;
        double c = xc*xc+ yc*yc + zc*zc- s.Radius()*s.Radius();
        double q1 = Quad1(a,b,c);
        double q2 = Quad2(a,b,c);
        
        // checks inside or if never hits
        if (q1 > 0) {
         i.t = q1;
         i.pos = R3Point(r.Start()+q1*r.Vector());
         i.norm = -1*R3Vector(i.pos.X()-s.Center().X(),i.pos.Y()-s.Center().Y(),i.pos.Z()-s.Center().Z());
            i.norm.Normalize();
        }
        else {
            i.t = q2;
            i.pos = R3Point(r.Start()+q2*r.Vector());
            i.norm = -1*R3Vector(i.pos.X()-s.Center().X(),i.pos.Y()-s.Center().Y(),i.pos.Z()-s.Center().Z());
            i.norm.Normalize();
        }
        
        return i;

    }
    
    
    // two hits
    else {
        i.intersected = true;
        
        //necessary dist
        double nec = sqrt(s.Radius()*s.Radius() - d*d);
        R3Vector v = R3Vector(s.Center()-r.Start());
        v.Project(r.Vector());
        double currt;
        if (r.Vector().X() != 0) currt = v.X()/r.Vector().X();
        else if (r.Vector().Y() != 0) currt = v.Y()/r.Vector().Y();
        else currt = v.Z()/r.Vector().Z();
        double distpert = r.Vector().Length();
        double realt = currt - nec/distpert;
        i.t = realt;
        i.pos = R3Point(r.Start()+realt*r.Vector());
        i.norm = R3Vector(i.pos.X()-s.Center().X(),i.pos.Y()-s.Center().Y(),i.pos.Z()-s.Center().Z());
        i.norm.Normalize();
    }
    return i;
}

// find intensity at a point given a light, using formulas from lecture
R3Rgb Intensity(R3Light light, R3Point p) {
    
    // directional light
    R3Rgb i0 = light.color;
    if (light.type == R3_DIRECTIONAL_LIGHT) {
        return i0;
    }
    double kc = light.constant_attenuation;
    double kl = light.linear_attenuation;
    double kq = light.quadratic_attenuation;
    double d = R3Distance(light.position, p);
    
    // point light
    if (light.type == R3_POINT_LIGHT) {

        return (i0 / (kc + kl*d + kq*d*d));
    }
    
    // spot light
    if (light.type == R3_SPOT_LIGHT) {
        R3Vector vectd = light.direction;
        R3Vector vectl = p - light.position;

        vectd.Normalize();
        vectl.Normalize();
        double dot = vectd.Dot(vectl);
        double acost = acos(dot);
     

        if (acost > light.angle_cutoff) return R3Rgb(0,0,0,0);
        R3Rgb ret = i0*(pow(dot,light.angle_attenuation))/(kc + kl*d + kq*d*d);
        return ret;
    }

    return i0;
}

// calculate diffuse intensity from lecture
R3Rgb DiffuseIntensity(R3Intersect sect, R3Light light, R3Material mat) {
    R3Vector l = (light.position-sect.pos);    
    if (light.type == R3_DIRECTIONAL_LIGHT) {
        l = -1*light.direction;
    }
    l.Normalize();
    sect.norm.Normalize();

    // if on wrong side of the material, gives 0 light
    if (sect.norm.Dot(l) < 0) return R2Pixel(0,0,0,0);
    R3Rgb total = mat.kd*Intensity(light,sect.pos)*(sect.norm.Dot(l));
    total.Clamp();

    return total;
}

// calcualte specular intensity from lecture
R3Rgb SpecularIntensity(R3Intersect sect, R3Light light, R3Material mat, R3Point camera) {

    R3Vector r = (sect.pos-light.position);
    r.Normalize();

    
    if (light.type == R3_DIRECTIONAL_LIGHT) {
        r = light.direction;
    }
    
    R3Plane p = R3Plane(sect.norm,0);
    
    r.Mirror(p);
    R3Vector v = camera - sect.pos;
    r.Normalize();
    v.Normalize();
    
    // if on wrong side, give 0 light
    if (v.Dot(r) < 0) return R2Pixel(0,0,0,0);
    R3Rgb total = mat.ks*pow((v.Dot(r)),mat.shininess)*Intensity(light,sect.pos);
    total.Clamp();

    return total;
}

// goes through scene, finds closest intersection and returns it
R3Intersect ComputeIntersect(R3Scene *scene, R3Node *node, R3Ray *r) {
    //Check for intersection with shape
    R3Intersect shape_intersection;
    
    
    // do transformation
    R3Matrix trans = node->transformation;
    trans.Invert();
    (*r).Transform(trans);
    
    // find intersection of current shape
    if (node->shape != NULL) {
        if (node->shape->type == R3_BOX_SHAPE) {
            shape_intersection = BoxIntersect(*r, *node->shape->box);
            shape_intersection.node = node;
        }
        else if (node->shape->type == R3_SPHERE_SHAPE) {
            shape_intersection = SphereInt(*r, *node->shape->sphere);
            shape_intersection.node = node;
            
        }
        else if (node->shape->type == R3_MESH_SHAPE) {
            shape_intersection = TriMeshIntersect(*r, *node->shape->mesh);
            shape_intersection.node = node;
        }
        else if (node->shape->type == R3_CYLINDER_SHAPE) {
            shape_intersection = CylinderIntersect(*r, *node->shape->cylinder);
            shape_intersection.node = node;
            
            
        }
        else if (node->shape->type == R3_CONE_SHAPE) {
            shape_intersection = ConeIntersect(*r, *node->shape->cone);
            shape_intersection.node = node;
        }
        
    }
    
    // if node was null, no closests intersection, go on to children and find closests
    else {
        R3Intersect closest = R3Intersect();
        closest.intersected = false;
        closest.t = std::numeric_limits<double>::infinity();
        
        for (unsigned int i = 0; i < node->children.size(); i++) {
            R3Intersect boxsect = BoxIntersect(*r,node->children[i]->bbox);
            if (boxsect.intersected && boxsect.t < closest.t) {
            R3Intersect childsect = ComputeIntersect(scene, node->children[i],r);
            if (childsect.intersected == true) {
                if (childsect.t < closest.t) {
                    
                    if (childsect.t > 1e-12) closest = childsect;
                }
            }
        }
        }
        
        // uninvert
        trans.Invert();
        closest.pos.Transform(trans);
        closest.norm.Transform(trans);
        (*r).Transform(trans);
        closest.t = R3Distance(closest.pos, (*r).Start());
        
        
        return closest;
    }
    
    // if no children, uninvert and return the shapes intersection
    if (node->children.size() == 0) {

        
        trans.Invert();
        (*r).Transform(trans);
        
        return shape_intersection;
    }
    
    // if it didnt intersect the bounding shape, no need to keep going
    if (shape_intersection.intersected == false) {
        trans.Invert();
        (*r).Transform(trans);
        
        return shape_intersection;
        
    }

    
    // if we are here, the node is a bounding shape for the other things, so the closest should be infinitely far away. Look through children and find closest if it exists.
    R3Intersect closest = R3Intersect();
    closest.intersected = false;
    closest.t = std::numeric_limits<double>::infinity();
    
    for (unsigned int i = 0; i < node->children.size(); i++) {
        R3Intersect boxsect = BoxIntersect(*r,node->children[i]->bbox);
        if (boxsect.intersected && boxsect.t < closest.t) {
        R3Intersect childsect = ComputeIntersect(scene, node->children[i],r);
        if (childsect.intersected == true) {
            if (childsect.t < closest.t) {
                
                if (childsect.t > 1e-12) closest = childsect;
            }
        }
        }
    }
    trans.Invert();
    closest.pos.Transform(trans);
    closest.norm.Transform(trans);
    (*r).Transform(trans);
    closest.t = R3Distance(closest.pos, (*r).Start());
    
    
    return closest;

}

// find the Vector the points towards a light given a point and a light source.
R3Vector LightVect(R3Light *light, R3Point p) {
    R3Vector v;
    if (light->type == R3_DIRECTIONAL_LIGHT) {
        v = light->direction;
        v.Normalize();
        return v;
    }
    v = p - light->position;
    v.Normalize();
    return v;
}

// finds a ray going towards a light given a point and light source
R3Ray LightRay(R3Light *light, R3Point p) {
    return R3Ray(p, -1*LightVect(light,p));
}

// finds the distance to a certain light
double LightDist(R3Light *light, R3Point p) {
    if  (light->type == R3_DIRECTIONAL_LIGHT) return std::numeric_limits<double>::infinity();
    return R3Distance(light->position,p);
}

// calcualte phong illumination as in lecture
R3Rgb Phong(R3Scene *scene, R3Ray r, R3Intersect hit) {
    R3Rgb total = R2Pixel(0,0,0,0);
    total += (hit.node->material->ka)*(scene->ambient);
    total += (hit.node->material->emission);
    for (int i = 0; i < scene->NLights(); i++) {
        R3Light *light = scene->Light(i);
        
        // calculate shadow, if shadow ray intersects and distance is closer than distance to light, cast shadow
        double shadown = 1.0;
        R3Ray shadowray = LightRay(light,hit.pos);
        R3Intersect shadowsect = ComputeIntersect(scene, scene->Root(), &shadowray);
        if (shadowsect.intersected == true) {
            if (shadowsect.t < LightDist(light,hit.pos))
            shadown = 0.0;

        }
        
        
        total += shadown*DiffuseIntensity(hit, *light, *hit.node->material);
        if (&hit != NULL) {
            
            if(light != NULL && hit.node->material != NULL) {
                total += shadown*SpecularIntensity(hit, *light, *hit.node->material,r.Start());
            }
            
            
        }
    }
    return total;

}

// computes the radiance given a ray and scene
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray ray, int max_depth, double currentindex) {
    //return R3Rgb(0,0,0,0);
    R3Intersect hit = ComputeIntersect(scene, scene->Root(),&ray);
    return ComputeRadiance(scene,ray,hit,max_depth, currentindex);
}

// computes radiance given a ray, scene and intersection, using phong illumination and other properties
R3Rgb ComputeRadiance(R3Scene *scene, R3Ray r, R3Intersect hit, int max_depth, double currentindex) {
    if (hit.intersected == false) {
     return scene->background;   
    }
    assert(hit.node != NULL);
    assert(hit.node->material != NULL);
    
    R3Rgb total = R2Pixel(0,0,0,0);
    total += Phong(scene, r, hit);
    
    if (max_depth > 0) {
        R3Rgb spec = ComputeRadiance(scene, SpecularRay(r,hit), max_depth-1,currentindex);
        total += (hit.node->material->ks)*spec;
        
        
        R3Rgb refract;
        if (currentindex == hit.node->material->indexofrefraction) {
            refract = ComputeRadiance(scene, RefractionRay(r,hit,currentindex,1), max_depth-1,1);
        }
        else {
        refract = ComputeRadiance(scene, RefractionRay(r,hit,currentindex,hit.node->material->indexofrefraction), max_depth-1,hit.node->material->indexofrefraction);
        }
        total += (hit.node->material->kt)*refract;
        
        R3Rgb trans = ComputeRadiance(scene, TransmissionRay(r,hit), max_depth-1,hit.node->material->indexofrefraction);
        total += (hit.node->material->kt)*trans;
        
    }
    return total;
}


// Makes picture by calling get color on reach ray, setting color of pixel.
R2Image *RenderImage(R3Scene *scene, int width, int height, int max_depth,
  int num_primary_rays_per_pixel, int num_distributed_rays_per_intersection)
{
  // Allocate  image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }
    
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            R3Ray r = ConstructRayThroughPixel(i,j,width,height,scene->Camera().eye,scene->Camera().towards,scene->Camera().up,scene->Camera().right,scene->Camera().xfov, scene->Camera().yfov);
            double currentindex = 1.0;
            R3Rgb radiance = ComputeRadiance(scene, r, max_depth, currentindex);
             image->SetPixel(i,j,radiance);
            
        }
    } 
  return image;
}


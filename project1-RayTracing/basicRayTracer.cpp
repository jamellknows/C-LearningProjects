#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#if defined __linux__ || defined __APPLE__
//"Compiled for Linux"
#else //Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8 
#endif
//This is a template for the class T 
/* Templates are features of C++
 that allows functions and classes to 
 operate with generic types, allowing a function or class 
 to work on many different data types */
 /*Template named T for the class Vec3 */ 

template<typename T>

class Vec3
{
    public:
    // variables x, y and z are of data type T
    T x, y,z;
      /* 
    Vec3 has 3 versions that it can be called
     returning different things depending on 
     the input parameters */
    Vec3(x, y, z) : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx, T yy, Tzz) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
  /* Vec3 accepts the address/refernce to a variable
  but only to normalize */
    Vec3& normalize()
    {   //Not too sure what is going on here but 
        // it is normalizing the vectors
        //I've never seen length2() used - not much documentation
        T nor2 = length2();
        if (nor2 > 0)
        {
            T invNor = 1/ sqrt(nor2);
            x *= invNor, y *= invNor, z *=invNor;
        }
        return *this; // pointer return of current object
    }
    // Various functions/methods performed to calculate products
    // relationships 
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator * (const Vec3<T> &v) const {return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    //&v is a refence variable that will be input 
    T dot(const Vec3<T> &v) const {return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T> operator += (const Vec3<T> &v){ x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T> operator *= (const Vec3<T> &v){ x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z;} // here length2() is defined
    T length() const {return sqrt(length2());}
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v) //writing to what is refernenced as os
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os; // writing the vectors of v to os
    }

};
/* typedef is used to create an additional name(alias)
for another data type, but does not create a new type */ 


/* Here we've created a new typedef/alias for Vec3 using floats
from template T */ 
typedef Vec3<float> Vec3f; 

//creating a class named sphere
class Sphere
{
    public:
        Vec3f center;
        float radius, radius2;
        Vec3f surfaceColour, emissionColour;
        float transparency, reflection;
        //Sphere function takes parameters, 
        // uses predefined functions 
        // returns nothing 
        //Computes a ray-sphere using the geometric solution
        Sphere (
            const Vec3f &c,
            const float &r,
            const Vec3f &sc,
            const float &refl = 0,
            const float &transp = 0,
            const Vec3f &ec = 0) :
            center(c), 
            radius(r), 
            radius2(r * r), 
            surfaceColour(sc), 
            emissionColour(ec),
            transparency(transp),
            reflection(refl)
            { /* empty */}
        
        //Calculates intersection
        bool intersect(const Vec3f &rayorig, 
                       const Vec3f &raydir,
                       float &t0, 
                       float &t1) const
                       {
                           Vec3f l = center - rayorig;
                           float tca = l.dot(raydir);
                           if(tca < 0) return false;
                           float thc = sqrt(radius2 - d2);
                           t0 = tca - thc;
                           t1 = tca + thc;

                           return true;
                       }
};

// This variable controls the maximum recursion depth 

#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

/*
This is the main function. It takes a ray an an argument 
(defined by it's origin and direction)
We test if this ray intersects any geometry in the scene. 
If the ray intersects an object, we compute the intersection
point, and shade this point using this information. 
Shading depends on the surface property (
is it transparent,reflective, diffuse). 
The function returns a colour for the ray, 
if the ray intersects an object that is the colour of the
object at the intersection point. 
Otherwise it returns a background colour.
*/

Vec3f trace(
    const Vec3f *rayorig, 
    const Vec3f &raydir, 
    const std::vector<Sphere> &spheres, 
    const int &depth)
{
        //if (raydir.length() != 1) std::cerr << "Error " <<raydir << std::endl;
        float tnear = INFINITY;
        const Sphere* sphere = NULL;

        //find intersection of this ray with the sphere in the scene
        for (unsigned i = 0; i < spheres.size(); ++i)
        {
            float t0 = INFINITY, t1 = INFINITY;
            if(spheres[i].intersect(rayorig, raydir, t0, t1))
            {
                if(t0 < 0) t0 = t1;
                if(t0 < tnear)
                {
                    tnear = t0;
                    sphere = &spheres[i];
                }
            }
        }
    

    //if there's no intersction return black or backgrounkd colour
    if (!sphere) return Vec3f(2);
    Vec3f surfaceColour = 0; //colour of the ray/surface of the object intersected by the ray 
    Vec3f phit = rayorig + raydir * tnear; //point of intersection
    Vec3f nhit = phit - sphere->center; // normal at the intersection point
    nhit.normalize(); // normalize normal direction
    /*If the normal and the view direction 
    are not opposite each other reverse the normal 
    direction. That also means we are inside the sphere 
    so set the inside bool to true. Finally reverse the sign
    of IdotN which we want positive
    */
   float bias = 1e-4; // add some bias to the point from which we will be tracing
   bool inside = false;
   if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
   if((sphere->transparency > 0)|| sphere->reflection > 0) && depth < MAX_RAY_DEPTH
   {
       float facingratio = - raydir.dot(nhit);
       //change the mix value to tweak the effect
       float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
       //compute reflection direction (no need to normalize vectors as already normalised )
       Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
       refldir.normalize();
       Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
       Vec3f refraction = 0;
       // if the sphere is also transparent compute refraction ray (transmission)
       if(sphere->transparency)
       {
           float ior = 1.1, eta = (inside) ? ior : 1/ ior; // are we inside or outside the surface?
           float cosi = -nhit.dot(raydir);
           float k = 1 - eta * eta * (1 - cosi * cosi);
           Vec3f refrdir = raydir * eta + nhit * (eta * cosi - sqrt(k));
           refrdir.normalsize();
           refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
       }
       //the result is a mix of reflection and refraction (If the sphere is transparent)
       surfaceColour = (
           reflection * fresneleffect +
           refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColour;

   }else 
   {
       //it's a diffuse object, no need to rayrtrace any further 
       for(unsigned i =0; i < spheres.size(); ++i)
       {
           if(spheres[i].emissionColour.x > 0)
           {
               //this is a light
               Vec3f transmission = 1;
               Vec3f lightDirection = spheres[i].center - phit;
               lightDirection.normalize();
               for(unsigned j = 0; j <spheres.size(); ++j)
               {
                   if(i != j )
                   {
                       float t0, t1;
                       if(shpheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1))
                       {
                           transmission = 0;
                           break;
                       }
                   }
               }

               surfaceColour += sphere->surfaceColour * transmission *
               std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColour;
           }

       }
   }

   return surfaceColour + sphere->emissionColour;
} 

/* 
Main rendering function. 
We compute a camera ray for each pixel of the image trace
it and return a colour. If the ray hits a sphere, we return 
the colour of the sphere at the intersection point, 
else we return the background colour.
*/
    
void render(const std::vector<Sphere> &spheres)
{
    unsigned width = 640, height = 480;
    // image size
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    //new image
    float invWidth = 1/ float(width), invHeight = 1/ float(height);
    // calculate normalized vectors
    float fov = 30, aspectratio = width / float(height);
    // field of view
    float angle = tan(M_PI * 0.5 * fov /180.);
    //Trace rays 
    for (unsigned y =0; y < height; ++y)
    {
        //casting 2d image
        float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
        float yy = (1 - 2* ((y + 0.5) * invHeight)) * angle;
        Vec3f raydir(xx, yy, -1);
        raydir.normalize();
        *pixel = trace(Vec3f(0), raydir, spheres, 0);
    }
//Save result to a PPm image (keep these flags if you compile under Windows)
std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
ofs << "P6\n" <<width << " " << height << "\n255\n";
for (unsigned i =0; i < width * height ; ++i)
{
    ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
    (unsigned char)(std::min(float(1), image[i].x) * 255) <<
    (unsigned char)(std::min(float(1), image[i].y) * 255) <<
    (unsigned char)(std::min(float(1), image[i].z) * 255);
}
ofs.close();
delete [] image;    
}

/* In the main function, 
we will create the scene which is composed
of 5 spheres and 1 light (which is also a sphere).
Then, once the scene description is complete
we render that scene, by calling the render() function.
*/

int main(int argc, char **argv)
{
    srand48(13);
    std::vector<Sphere> spheres;
    //position, radius, surface colour, reflectivity, transparency, emission colour
    spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, 0, -20), 4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vec3f( 5.0, -1, -15), 2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( 5,0, 0, -25), 3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( -5.5, 0, -15), 3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
    //light
    spheres.push_back(Sphere(Vec3f( 0.0, 20, 30), 3, Vec3f(0.00, 0.00, 0.0), 1, 0.0));
    render(spheres);

    return 0;

}
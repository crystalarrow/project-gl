/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

vec3 color;
struct Vertex
{
	vec4 position;
	vec3 color;
};
vector<Vertex> verticies;
Vertex vertex;
MGLpoly_mode polymode;
struct Triangle
{
	Vertex A, B, C;
};
vector<Triangle> triangles;
vector<mat4> stacks[2] = {{mat4()},{mat4()}};
MGLmatrix_mode matrixmode;
MGLfloat* zbuffer = 0;



float Calculate_Area(float ax, float ay, float bx, float by, float cx, float cy)
{
	return ax * (by - cy) + ay * (cx - bx) + (bx * cy - by * cx);
}

bool clipped(float x, float y, float z)
{
	return (x < -1 || x > 1 || y < -1 || y > 1 || z < -1 || z > 1);
}

void Rasterize_Triangle(const Triangle& tri, 
						int width, 
						int height, 
						MGLpixel* data)
{
	// point A
	float ai = ((1 + tri.A.position[0]) * width) / 2 - 0.5;
	float aj = ((1 + tri.A.position[1]) * height) / 2 - 0.5;
	// point B
	float bi = ((1 + tri.B.position[0]) * width) / 2 - 0.5;
	float bj = ((1 + tri.B.position[1]) * height) / 2 - 0.5;
	// point C
	float ci = ((1 + tri.C.position[0]) * width) / 2 - 0.5;
	float cj = ((1 + tri.C.position[1]) * height) / 2 - 0.5;

	// calculates area for entire triangle
	float abc_area = Calculate_Area(ai, aj, bi, bj, ci, cj);

	for(int i = max(min(min((int)ai, (int)bi), (int)ci), 0); i <= min(max(max((int)ai,(int) bi), (int)ci), width-1); i++)
	{	
		for(int j = max(min(min((int)aj, (int)bj),(int) cj), 0); j <= min(max(max((int)aj, (int)bj), (int)cj), height-1); j++)
		{
			vec3 p(i,j,0); // test point

			/* convert to barycentric */
			// pbc
			float alpha = Calculate_Area(p[0], p[1], bi, bj, ci, cj) / abc_area;
			// apc
			float beta = Calculate_Area(ai, aj, p[0], p[1], ci, cj) / abc_area;
			// abp
			float gamma = Calculate_Area(ai, aj, bi, bj, p[0], p[1]) / abc_area;

			// checks if the test point is inside the triangle
			if(alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				float x = 2.0 * i / width - 1.0;
				float y = 2.0 * j / height - 1.0;
				float z = alpha*tri.A.position[2] + beta*tri.B.position[2] + gamma*tri.C.position[2]; 
				if(!clipped(x,y,z) && z < zbuffer[i + j*width])
				{
					zbuffer[i + j*width] = z;
					
					vec3 pcolor = tri.A.color*alpha + tri.B.color*beta + tri.C.color*gamma;
					data[i + j*width] = Make_Pixel(pcolor[0]*255, pcolor[1]*255, pcolor[2]*	255);
				}
			}
		}
	}
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize widthmin,
                   MGLsize height,
                   MGLpixel *data)
{
	if (zbuffer != 0) free(zbuffer);
	zbuffer = (MGLfloat*) malloc(widthmin * height * sizeof(MGLfloat));


	for(unsigned int i = 0; i < widthmin*height; i++)
	{	
   		data[i] = Make_Pixel(0,0,0);
   		zbuffer[i] = 2;
	}

	for(unsigned int i = 0; i < triangles.size(); i++)
	{
		Rasterize_Triangle(triangles[i], widthmin, height, data);
	}

	triangles.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	polymode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(polymode == MGL_TRIANGLES)
	{
		for(unsigned int i = 0; i < verticies.size(); i+=3)
		{
			Triangle t ={.A = verticies[i], .B = verticies[i + 1], .C = verticies[i + 2]};
			triangles.push_back(t);

		}
	}
	if(polymode == MGL_QUADS)
	{
		for(unsigned int i = 0; i < verticies.size() ; i+=4)
		{
			Triangle t1 ={.A = verticies[i], .B = verticies[i + 1], .C = verticies[i + 2]};
			Triangle t2 ={.A = verticies[i], .B = verticies[i + 2], .C = verticies[i + 3]};
			triangles.push_back(t1);
			triangles.push_back(t2);
		}	
	}

	verticies.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{

	mglVertex3(x, y, 0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{ 
	vec4 p(x, y, z, 1);
	
	p = stacks[MGL_MODELVIEW].back() * p; 
	p = stacks[MGL_PROJECTION].back() * p; 
	p /= p[3];
	verticies.push_back({p, color});
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	matrixmode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	stacks[matrixmode].push_back(stacks[matrixmode].back());
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	stacks[matrixmode].pop_back();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	mat4 identity;

	identity.make_zero();
	identity(0,0) = 1;
	identity(1,1) = 1;
	identity(2,2) = 1;
	identity(3,3) = 1;

	stacks[matrixmode].back() = identity;
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    for(int i = 0; i < 16; i++)
        stacks[matrixmode].back().values[i] = matrix[i];
}
/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	mat4 m;

    for(int i = 0; i < 16; i++)
        m.values[i] = matrix[i];	

    stacks[matrixmode].back() = stacks[matrixmode].back() * m;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat4 translate;
	translate.make_zero();

	translate(0,0) = 1;
	translate(1,1) = 1;
	translate(2,2) = 1;
	translate(3,3) = 1;
	translate(0,3) = x;
	translate(1,3) = y;
	translate(2,3) = z;

	mglMultMatrix(translate.values);
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
	mat4 rotate;
	rotate.make_zero();

	// need to convert angle into degrees
	float c = cos(angle * M_PI / 180);
	float s = sin(angle * M_PI / 180);

	// normalize the values
	float magnitude = sqrt(x*x + y*y + z*z);
	float ax = x / magnitude;
	float ay = y / magnitude;
	float az = z / magnitude;

	rotate(0,0) = ax*ax*(1-c)+c;
	rotate(0,1) = ax*ay*(1-c)-az*s;
	rotate(0,2) = ax*ay*(1-c)+ay*s;
	rotate(1,0) = ay*ax*(1-c)+az*s;
	rotate(1,1) = ay*ay*(1-c)+c;
	rotate(1,2) = ay*az*(1-c)-ax*s;
	rotate(2,0) = ax*az*(1-c)-ay*s;
	rotate(2,1) = ay*az*(1-c)+ax*s;
	rotate(2,2) = az*az*(1-c)+c;
	rotate(3,3) = 1;

	mglMultMatrix(rotate.values);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 scale;
	scale.make_zero();

	scale(0,0) = x;
	scale(1,1) = y;
	scale(2,2) = z;
	scale(3,3) = 1;

	mglMultMatrix(scale.values);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	mat4 frustum;
	frustum.make_zero();
	frustum(0,0) = 2 * near / (right - left);
	frustum(1,1) = 2 * near/ (top - bottom);
	frustum(0,2) = (right + left) / (right - left);
	frustum(1,2) = (top + bottom) / (top - bottom);
	frustum(2,2) = -(far + near) / (far - near);
	frustum(2,3) = -(2 * far * near) / (far - near);
	frustum(3,2) = -1;

	mglMultMatrix(frustum.values);
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	mat4 ortho;
	ortho.make_zero();
	ortho(0,0) = 2 / (right - left);
	ortho(1,1) = 2 / (top - bottom);
	ortho(2,2) = -2 / (far - near);
	ortho(0,3) = -(right + left) / (right - left); 
	ortho(1,3) = -(top + bottom) / (top - bottom);
	ortho(2,3) = -(far + near) / (far - near);
	ortho(3,3) = 1;

	mglMultMatrix(ortho.values);
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	color = vec3(red, green, blue);
}

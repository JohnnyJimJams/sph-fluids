/* $Id: Vector.cpp,v 1.1.1.1 2005/05/01 15:18:55 ovidiom Exp $ */

#include "Vector.h"


Vector3f operator +(const Vector3f &v, float s)
{
	return (Vector3f(v.x + s, v.y + s, v.z + s));
}

Vector3f operator +(float s, const Vector3f &v)
{
	return (Vector3f(s + v.x, s + v.y, s + v.z));
}

Vector3f operator +(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x + v.x, u.y + v.y, u.z + v.z));
}

Vector3f operator -(const Vector3f &v, float s)
{
	return (Vector3f(v.x - s, v.y - s, v.z - s));
}

Vector3f operator -(float s, const Vector3f &v)
{
	return (Vector3f(s - v.x, s - v.y, s - v.z));
}

Vector3f operator -(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x - v.x, u.y - v.y, u.z - v.z));
}

Vector3f operator *(const Vector3f &v, float s)
{
	return (Vector3f(v.x * s, v.y * s, v.z * s));
}

Vector3f operator *(float s, const Vector3f &v)
{
	return (Vector3f(s * v.x, s * v.y, s * v.z));
}

Vector3f operator *(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x * v.x, u.y * v.y, u.z * v.z));
}

Vector3f operator /(const Vector3f &v, float s)
{
	float	inv = 1.0f / s;

	return (Vector3f(v.x * inv, v.y * inv, v.z * inv));
}

Vector3f operator /(float s, const Vector3f &v)
{
	return (Vector3f(s / v.x, s / v.y, s / v.z));
}

Vector3f operator /(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x / v.x, u.y / v.y, u.z / v.z));
}

Vector3f operator -(const Vector3f &v)
{
	return (Vector3f(-v.x, -v.y, -v.z));
}

float dot(const Vector3f &u, const Vector3f &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

Vector3f cross(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.y * v.z - v.y * u.z,
	                 u.z * v.x - u.x * v.z,
	                 u.x * v.y - u.y * v.x));
}

float length(const Vector3f &v)
{
	return ((float)sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

Vector3f normalize(const Vector3f &v)
{
	return (v / length(v));
}


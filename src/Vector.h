/* $Id: Vector.h,v 1.1.1.1 2005/05/01 15:18:55 ovidiom Exp $ */

#ifndef _VECTOR3F_H_
#define _VECTOR3F_H_

#include <cmath>


struct Vector3f
{
	float	x, y, z;

	inline Vector3f()
	{
		x = y = z = 0.0f;
	}
	
	inline Vector3f(float x, float y, float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	
	inline Vector3f(float xyz)
	{
		x = y = z = xyz;
	}

	inline Vector3f(const float *xyzArr)
	{
		x = xyzArr[0];
		y = xyzArr[1];
		z = xyzArr[2];
	}
	
	inline operator const float *() const
	{
		return ((const float *)&x);
	}
	
	inline float &operator[](unsigned int idx)
	{
		return (*(((float *)&x) + idx));
	}
	
	inline void operator +=(float s)
	{
		x += s;
		y += s;
		z += s;
	}
	
	inline void operator +=(const Vector3f &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}
	
	inline void operator -=(float s)
	{
		x -= s;
		y -= s;
		z -= s;
	}
	
	inline void operator -=(const Vector3f &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}
	
	inline void operator *=(float s)
	{
		x *= s;
		y *= s;
		z *= s;
	}
	
	inline void operator *=(const Vector3f &v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
	}
	
	inline void operator /=(float s)
	{
		float	inv = 1.0f / s;

		x *= inv;
		y *= inv;
		z *= inv;
	}

	inline void operator /=(const Vector3f &v)
	{
		x /= v.x;
		y /= v.y;
		z /= v.z;
	}
};


inline Vector3f operator +(const Vector3f &v, float s)
{
	return (Vector3f(v.x + s, v.y + s, v.z + s));
}

inline Vector3f operator +(float s, const Vector3f &v)
{
	return (Vector3f(s + v.x, s + v.y, s + v.z));
}

inline Vector3f operator +(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x + v.x, u.y + v.y, u.z + v.z));
}

inline Vector3f operator -(const Vector3f &v, float s)
{
	return (Vector3f(v.x - s, v.y - s, v.z - s));
}

inline Vector3f operator -(float s, const Vector3f &v)
{
	return (Vector3f(s - v.x, s - v.y, s - v.z));
}

inline Vector3f operator -(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x - v.x, u.y - v.y, u.z - v.z));
}

inline Vector3f operator *(const Vector3f &v, float s)
{
	return (Vector3f(v.x * s, v.y * s, v.z * s));
}

inline Vector3f operator *(float s, const Vector3f &v)
{
	return (Vector3f(s * v.x, s * v.y, s * v.z));
}

inline Vector3f operator *(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x * v.x, u.y * v.y, u.z * v.z));
}

inline Vector3f operator /(const Vector3f &v, float s)
{
	float	inv = 1.0f / s;

	return (Vector3f(v.x * inv, v.y * inv, v.z * inv));
}

inline Vector3f operator /(float s, const Vector3f &v)
{
	return (Vector3f(s / v.x, s / v.y, s / v.z));
}

inline Vector3f operator /(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.x / v.x, u.y / v.y, u.z / v.z));
}

inline Vector3f operator -(const Vector3f &v)
{
	return (Vector3f(-v.x, -v.y, -v.z));
}

inline float dot(const Vector3f &u, const Vector3f &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

inline Vector3f cross(const Vector3f &u, const Vector3f &v)
{
	return (Vector3f(u.y * v.z - v.y * u.z,
	                 u.z * v.x - u.x * v.z,
	                 u.x * v.y - u.y * v.x));
}

inline float length(const Vector3f &v)
{
	return ((float)sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

inline Vector3f normalize(const Vector3f &v)
{
	return (v / length(v));
}



#endif /* _VECTOR3F_H_ */


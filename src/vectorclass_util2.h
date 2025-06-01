#ifndef VECTORCLASS_UTIL2_H
#define VECTORCLASS_UTIL2_H  

typedef Vec4fb Vec2fb;

#if 1
typedef Vec4ui Vec2ui;
#else
class Vec2ui {
protected:
    __m128i xmm; // Integer vector

public:
    // Default constructor:
    Vec2ui() {
    }
    // Constructor to broadcast the same value into all elements:
    Vec2ui(uint32_t i) {
        xmm = _mm_setr_epi32((int32_t)i,(int32_t)i,0,0);
    }
    // Constructor to build from all elements:
    Vec2ui(uint32_t i0, uint32_t i1) {
        xmm = _mm_setr_epi32((int32_t)i0, (int32_t)i1, 0, 0);
    }
    // Constructor to convert from type __m128i used in intrinsics:
    Vec2ui(__m128i const & x) {        
        xmm = x;
	cutoff2;
    }
    // Assignment operator to convert from type __m128i used in intrinsics:
    Vec2ui & operator = (__m128i const & x) {
        xmm = x;
	cutoff2();
        return *this;
    }
    // Member function to load from array (unaligned)
    Vec2ui & load(void const * p) {
        return load_partial2(p);
    }
    // Member function to load from array (aligned)
    Vec2ui & load_a(void const * p) {
        return load_partial2(p);
    }
    // Member function extract a single element from vector
    uint32_t extract(int index) const {
        uint32_t i[4];
	_mm_storeu_epi32(i,xmm);
	return i[index & 3];
    }
    // Extract a single element. Use store function if extracting more than one element.
    // Operator [] can only read an element, not write.
    uint32_t operator [] (int index) const {
        return extract(index);
    }
};
#endif

#if 0
typedef Vec4f Vec2f;
#else
class Vec2f {
protected:
    __m128 xmm; // Float vector
public:
    // Default constructor:
    Vec2f() {
    }
    // Constructor to broadcast the same value into all elements:
    Vec2f(float f) {
        xmm = _mm_setr_ps(f,f,0.0f,0.0f); 
    }
    // Constructor to build from all elements:
    Vec2f(float f0, float f1) {
        xmm = _mm_setr_ps(f0, f1, 0.f, 0.f); 
    }
    // Constructor to convert from type __m128 used in intrinsics:
    Vec2f(__m128 const & x) {
        xmm = x;
	cutoff2();
    }
    // Assignment operator to convert from type __m128 used in intrinsics:
    Vec2f & operator = (__m128 const & x) {
        xmm = x;
        cutoff2();	
        return *this;
    }
    // Type cast operator to convert to Vec4f
    operator Vec4f() const {
        return xmm;
    }
    // Type cast operator to convert to __m128 used in intrinsics
    operator __m128() const {
        return xmm;
    }
    // Member function to load from array (unaligned)
    Vec2f & load(float const * p) {
        return load_partial2(p);
    }

    Vec2f & load_a(float const * p) {
        return load_partial2(p);       
    }

    // Member function to store into array (unaligned)
    void store(float * p) const {
        store_partial2(p);
    }

    void store_a(float * p) const {
        store_partial2(p);
    }

    // Partial load. Load n elements and set the rest to 0
    Vec2f & load_partial2(float const * p) {
        xmm = _mm_castpd_ps(_mm_load_sd((double const*)p));
        return *this;
    }

    // Partial store. Store n elements
    void store_partial2(float * p) const {
            _mm_store_sd((double*)p, _mm_castps_pd(xmm)) ;
    }

    Vec2f & cutoff2(void) {
        static const union {        
            int32_t i[8];
            float   f[8];
        } mask = {{1,-1,-1,-1,0,0,0,0}};
        xmm = _mm_and_ps(xmm, Vec4f().load(mask.f + 2));
        return *this;
    }

    // Member function extract a single element from vector
    float extract(int index) const {
        float x[4];
	_mm_storeu_ps(x, xmm);
        return x[index & 3];
    }
    // Extract a single element. Use store function if extracting more than one element.
    // Operator [] can only read an element, not write.
    float operator [] (int index) const {
        return extract(index);
    }
};

static inline Vec2f operator + (Vec2f const & a, Vec2f const & b) {
    return _mm_add_ps(a, b);
}
static inline Vec2f operator + (Vec2f const & a, float b) {
    return a + Vec2f(b);
}
static inline Vec4f operator + (float a, Vec2f const & b) {
    return Vec2f(a) + b;
}

static inline Vec2f operator - (Vec2f const & a, Vec2f const & b) {
    return _mm_sub_ps(a, b);
}

static inline Vec2f operator - (float a, Vec2f const & b) {
    return Vec2f(a) - b;
}

static inline Vec2f operator * (Vec2f const & a, Vec2f const & b) {
    return _mm_mul_ps(a, b);
}

static inline Vec2f & operator *= (Vec2f & a, Vec2f const & b) {
    a = a * b;
    return a;
}

static inline Vec2f mul_add(Vec2f const & a, Vec2f const & b, Vec2f const & c) {
    return a * b + c;
}

static inline Vec2f operator / (Vec2f const & a, Vec2f const & b) {
    Vec4f a_=(Vec4f)a,
	  b_=(Vec4f)b + Vec4f(0,0,1,1);
    Vec4f c_ = a_ / b_;
    return (Vec2f)c_;
}
static inline Vec2f operator / (float a, Vec2f const & b) {
    return Vec2f(a) / b;
}

static inline Vec2f abs(Vec2f const & a) {
    __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
    return _mm_and_ps(a,mask);
}

static inline Vec2f min(Vec2f const & a, Vec2f const & b) {
    return _mm_min_ps(a,b);
}

inline void vload_from_uint16_a( Vec2f &z, const uint16_t *p )
{
  Vec8us y(p[0],p[1],0,0,0,0,0,0);
  Vec4ui w = _mm_unpacklo_epi16(y, _mm_set1_epi16(0));
  z = _mm_cvtepi32_ps(w);
}

inline void vstream(float * p, Vec2f x) 
{
  x.store(p);
}

#endif

#endif  // VECTORCLASS_UTIL2_H


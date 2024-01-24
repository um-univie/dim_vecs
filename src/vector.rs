use rand::Rng;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Sub};
use num_traits::real::Real;

/// A 3D vector represented by its x, y, and z components.
///
/// # Example
///
/// ```
/// use crate::dimensional_vectors::vector::Vector3D;
///
/// let v1 = Vector3D::new(1.0, 2.0, 3.0);
/// let v2 = Vector3D::new(4.0, 5.0, 6.0);
/// let v3 = v1 + v2;
/// let v4 = v1 * 2.0;
/// let v5 = v1 / 2.0;
/// println!("Vector sum: {:?}", v3);
/// println!("Vector scalar multiplication: {:?}", v3);
/// println!("Vector scalar division: {:?}", v3);
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Vector3D<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl <T> Add<Vector3D<T>> for Vector3D<T> where T: Real {
    type Output = Vector3D<T>;

    fn add(self, other: Vector3D<T>) -> Vector3D<T> {
        Vector3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl <T> Mul<T> for Vector3D<T> where T: Real {
    type Output = Vector3D<T>;

    fn mul(self, scalar: T) -> Vector3D<T> {
        Vector3D {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

impl <T> Div<T> for Vector3D<T> where T: Real  {
    type Output = Vector3D<T>;

    fn div(self, scalar: T) -> Self::Output {
        Vector3D {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }
}
impl <T> Sub<Vector3D<T>> for Vector3D<T> where T: Real + From<f64> + Into<f64> {
    type Output = Vector3D<T>;
    fn sub(self, rhs: Vector3D<T>) -> Self::Output {
        Vector3D::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl <T> AddAssign for Vector3D<T> where T: Real + AddAssign {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

// Implement the AddAssign trait for Vector
impl <T> MulAssign<T> for Vector3D<T> where T: Real + MulAssign {
    fn mul_assign(&mut self, scalar: T) {
        self.x = self.x * scalar;
        self.y = self.y * scalar;
        self.z = self.z * scalar;
    }
}

impl <T> Default for Vector3D<T> where T: Real + From<f64> + Into<f64> {
    /// Creates a default `Vector` being the zero vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v: Vector3D<f64> = Vector3D::default();
    /// println!("Vector: {:?}", v);
    /// assert_eq!(v.length(),0.0)
    /// ```
    fn default() -> Vector3D<T> {
        Vector3D::new(0.0.into(), 0.0.into(), 0.0.into())
    }
}

impl <T> Vector3D<T> where T: Real + From<f64> + Into<f64> {
    /// Creates a `Vector` from an array of 3 elements, representing the x, y, and z components.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let array = [1.0, 2.0, 3.0];
    /// let v = Vector3D::from_vec(&array);
    /// println!("Vector: {:?}", v);
    /// ```
    pub fn from_vec(array: &[T; 3]) -> Vector3D<T> {
        Vector3D {
            x: array[0],
            y: array[1],
            z: array[2],
        }
    }
    /// Creates a `Vector` from x, y, and z components.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v = Vector3D::new(1.0, 2.0, 3.0);
    /// println!("Vector: {:?}", v);
    /// ```
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }
    pub fn x() -> Self {
        Self::new(1.0.into(), 0.0.into(), 0.0.into())
    }
    pub fn y() -> Self {
        Self::new(0.0.into(), 1.0.into(), 0.0.into())
    }
    pub fn z() -> Self {
        Self::new(0.0.into(), 0.0.into(), 1.0.into())
    }
    /// Calculates the difference between two `Vector` objects, this is equivalent to the vector
    /// pointing from `other` to `self`.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v1 = Vector3D::new(1.0, 1.0, 1.0);
    /// let v2 = Vector3D::new(1.0, 1.0, -1.0);
    ///
    /// let v_diff = v1.difference(&v2);
    /// assert_eq!(v_diff, Vector3D::new(0.0, 0.0, 2.0));
    /// assert_eq!(v_diff.length(), 2.0);
    /// ```
    pub fn difference(&self, other: &Self) -> Self {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        Self::new(dx, dy, dz)
    }
    /// Calculates the length of a vector.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    ///     let v = Vector3D::new(1.0, 1.0, 1.0);
    ///     println!("Vector: {:?}", v.length());
    ///     assert_eq!(v.length(), 3.0_f64.sqrt())
    ///
    /// ```
    pub fn length(&self) -> f64 {
        (self.x.into().powi(2) + self.y.into().powi(2) + self.z.into().powi(2)).sqrt()
    }
    /// Calculates the dot product of two vectors.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v1 = Vector3D::new(1.0, 2.0, 3.0);
    /// let v2 = Vector3D::new(4.0, 5.0, 6.0);
    /// println!("Vector: {:?}", v1.dot(&v2));
    ///
    /// ```
    pub fn dot(&self, other: &Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    /// Calculates the distance between to vectors
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v1 = Vector3D::new(1.0, 2.0, 3.0);
    /// let v2 = Vector3D::new(4.0, 5.0, 6.0);
    /// println!("Vector: {:?}", v1.distance(&v2));
    ///
    /// ```
    pub fn distance(&self, other: &Self) -> f64 {
        self.difference(other).length()
    }
    /// Returns the squared distance between two vectors
    ///
    /// # Example
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    /// let v1 = Vector3D::new(1.0, 2.0, 3.0);
    /// let v2 = Vector3D::new(4.0, 5.0, 6.0);
    /// assert_eq!(v1.distance_squared(&v2),27.0);
    /// ```
    pub fn distance_squared(&self, other: &Self) -> T {
        let difference = *other - *self;
        difference.x.powi(2) + difference.y.powi(2) + difference.z.powi(2)
    }
    /// Generates a random unit vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let v: Vector3D<f64> = Vector3D::random_unit_vector();
    /// println!("Vector: {:?}", v);
    /// assert_eq!(v.length().ceil(),1.0)
    /// ```
    ///
    pub fn random_unit_vector() -> Self {
        let mut rng = rand::thread_rng();
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let phi = rng.gen_range(0.0..2.0 * std::f64::consts::PI);

        let x = theta.sin() * phi.cos();
        let y = theta.sin() * phi.sin();
        let z = theta.cos();

        Vector3D::new(x.into(), y.into(), z.into())
    }
    /// Calculates the cross product of two vectors and returns a new vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let vec1: Vector3D<f64> = Vector3D::x();
    /// let vec2: Vector3D<f64> = Vector3D::y();
    /// assert_eq!(vec1.cross(&vec2),Vector3D::z())
    /// ```
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    /// Returns the angle between two vectors in radians.
    ///
    /// # Examples
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let vec_a = Vector3D::new(1.0, 0.0, 0.0);
    /// let vec_b = Vector3D::new(0.0, 1.0, 0.0);
    ///
    /// let angle = vec_a.angle_between(&vec_b);
    /// assert_eq!(angle, std::f64::consts::FRAC_PI_2);
    /// ```
    pub fn angle_between(&self, other: &Vector3D<T>) -> f64 {
        let lengths_product = self.length() * other.length();
        if lengths_product == 0.0 {
            0.0
        } else {
            let angle_cosine: f64 = (self.dot(other).into() / (lengths_product)).clamp(-1.0, 1.0);
            angle_cosine.acos()
        }
    }
    /// Calculates a normal vector for a given input vector
    ///
    /// # Panics
    /// This function does not panic as the case of a zero length vector is handeled.
    ///
    /// # Example
    ///
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let vec1 = Vector3D::new(1.0,2.0,3.0);
    /// assert_eq!(vec1.normalize().length().ceil(),1.0)
    /// ```
    pub fn normalize(&self) -> Self {
        if self.length() == 0.0 {
            *self
        } else {
            // TODO: This is a bit of a hack, but it works for now
            *self / self.length().into()
        }
    }
    /// Calculates the angle given three points in the order, start, middle, end
    ///
    /// # Panics
    /// Function does not panic
    ///
    /// # Example
    /// ```
    /// use crate::dimensional_vectors::vector::Vector3D;
    ///
    /// let start_point = Vector3D::new(1.0,0.0,0.0);
    /// let middle_point = Vector3D::new(0.0,0.0,0.0);
    /// let end_point = Vector3D::new(0.0,1.0,0.0);
    /// assert_eq!(start_point.angle_between_points(&middle_point,&end_point),1.5707963267948966);
    /// ```
    pub fn angle_between_points(&self, middle_point: &Vector3D<T>, end_point: &Vector3D<T>) -> f64 {
        let v1 = *middle_point - *self;
        let v2 = *end_point - *middle_point;
        v1.angle_between(&v2)
    }
    pub fn as_tuple(&self) -> (T, T, T) {
        (self.x, self.y, self.z)
    }
    pub fn as_array(&self) -> [T; 3] {
        [self.x, self.y, self.z]
    }
    pub fn as_vec(&self) -> Vec<T> {
        vec![self.x, self.y, self.z]
    }
}


use crate::vector::Vector3D;
#[test]
fn test_vector() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);
    let v2 = Vector3D::new(1.0, 2.0, 3.0);
    assert_eq!(v1, v2);
}
#[test]
fn test_vector_add() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);
    let v2 = Vector3D::new(1.0, 2.0, 3.0);
    let v3 = v1 + v2;
    assert_eq!(v3, Vector3D::new(2.0, 4.0, 6.0));
}
#[test]
fn test_vector_sub() {
    let v1 = Vector3D::new(1.0, 2.0, 3.0);
    let v2 = Vector3D::new(1.0, 2.0, 3.0);
    let v3 = v1 - v2;
    assert_eq!(v3, Vector3D::new(0.0, 0.0, 0.0));
}

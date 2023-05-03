
import logging
import numpy as np



def get_perp(v):
    """Return a vector perpendicular to the given vector.
    Will return the vector [0, 1, 0] if v is parallel to [0, 0, 1] or v is the zero vector.
    
    Args:
        v (np.ndarray): 3D vector.
        
    Returns:
        np.ndarray: 3D vector perpendicular to the given vector.
            
    """
    if np.allclose(v, 0):
        return np.array([0, 1, 0])
    v = v / np.linalg.norm(v)
    if np.allclose(np.abs(v), [0, 0, 1]):
        return np.array([0, 1, 0])
    cross = np.cross(v, [0, 0, 1])
    return cross / np.linalg.norm(cross)



def rigid_transform(A, B, sanity_check_tolerance=None):
    """Find the rigid transformation that aligns A to B.

    Args:
        A (np.ndarray): 3xN matrix of N 3D points.
        B (np.ndarray): 3xN matrix of N 3D points.

    Returns:
        tuple(np.ndarray, np.ndarray): Rotation matrix, translation vector.
    
    """

    # taken from https://github.com/nghiaho12/rigid_transform_3D/blob/master/rigid_transform_3D.py
    assert A.shape == B.shape

    num_rows, num_cols = A.shape
    if num_cols != 3:
        raise Exception(f"matrix A is not Nx3, it is {num_rows}x{num_cols}")

    num_rows, num_cols = B.shape
    if num_cols != 3:
        raise Exception(f"matrix B is not Nx3, it is {num_rows}x{num_cols}")

    A = A.T
    B = B.T

    if sanity_check_tolerance is not None:
        pairwise_distances_1 = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
        pairwise_distances_2 = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
        if np.max(np.abs(pairwise_distances_1 - pairwise_distances_2).flatten()) > sanity_check_tolerance:
            print(f"Pairwise distances deviate by {np.max(np.abs(pairwise_distances_1 - pairwise_distances_2).flatten())}. It is likely the points were not chosen well.")


    # find mean column wise
    centroid_A = np.mean(A, axis=1)
    centroid_B = np.mean(B, axis=1)

    # ensure centroids are 3x1
    centroid_A = centroid_A.reshape(-1, 1)
    centroid_B = centroid_B.reshape(-1, 1)

    # subtract mean
    Am = A - centroid_A
    Bm = B - centroid_B

    H = Am @ np.transpose(Bm)

    # find rotation
    U, S, Vt = np.linalg.svd(H)
    u_det = np.linalg.det(U)
    v_det = np.linalg.det(Vt)
    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        logging.debug("det(R) < R, reflection detected!, correcting for it ...")
        Vt[2,:] *= -1
        R = Vt.T @ np.diag([1, 1, u_det * v_det]) @ U.T

    t = -R @ centroid_A + centroid_B

    return R, t.flatten()
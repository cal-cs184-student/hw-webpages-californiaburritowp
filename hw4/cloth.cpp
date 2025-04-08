#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
// Create grid of point masses and store them in row-major order.
    for (int i = 0; i < num_height_points; i++) {
        // j is our x axis index (column #), i is our y axis index (row #)
        for (int j = 0; j < num_width_points; j++) {
            // Calculate normalized coordinates 
            double u = static_cast<double>(j) / (num_width_points - 1);
            double v = static_cast<double>(i) / (num_height_points - 1);

            Vector3D pos;
            // Check orientation: HORIZONTAL cloth lies in x-z plane with a fixed y value.
            if (orientation == HORIZONTAL) {
                double x = u * width;
                double y = 1.0;          // Fixed y coordinate for horizontal cloth.
                double z = v * height;
                pos = Vector3D(x, y, z);
            }
            else { // VERTICAL cloth lies in x-y plane, with a small random offset for z.
                double x = u * width;
                double y = v * height;
                // Generate a small random offset for z between -1/1000 and 1/1000.
                // rand() / RAND_MAX is a value in the range 0 and 1. 
                double rand_offset = ((rand() / static_cast<double>(RAND_MAX)) * (2.0 / 1000.0))
                    - (1.0 / 1000.0);
                double z = rand_offset;
                pos = Vector3D(x, y, z);
            }

            // Check if the mass should be pinned based on the pinned vector.
            bool isPinned = false;
            if (i < pinned.size()) {
                for (int pinned_index : pinned[i]) {
                    if (pinned_index == j) {
                        isPinned = true;
                        break;
                    }
                }
            }
            
            // Add the new point mass to the end of the vector
            point_masses.emplace_back(pos, isPinned);
        }
    }

    // Now that all point masses are created, create springs to connect them.
    for (int i = 0; i < num_height_points; i++) {
        for (int j = 0; j < num_width_points; j++) {
            // Compute the 1D index from the 2D grid coordinates.
            int idx = i * num_width_points + j;

            // STRUCTURAL springs: left and above.
            if (j > 0) {
                int left_idx = i * num_width_points + (j - 1);
                springs.emplace_back(&point_masses[idx], &point_masses[left_idx], STRUCTURAL);
            }
            if (i > 0) {
                int above_idx = (i - 1) * num_width_points + j;
                springs.emplace_back(&point_masses[idx], &point_masses[above_idx], STRUCTURAL);
            }

            // SHEARING springs: upper-left and upper-right.
            if (i > 0 && j > 0) {
                int diagLeft_idx = (i - 1) * num_width_points + (j - 1);
                springs.emplace_back(&point_masses[idx], &point_masses[diagLeft_idx], SHEARING);
            }
            if (i > 0 && j < num_width_points - 1) {
                int diagRight_idx = (i - 1) * num_width_points + (j + 1);
                springs.emplace_back(&point_masses[idx], &point_masses[diagRight_idx], SHEARING);
            }

            // BENDING springs: two masses away horizontally and vertically.
            if (j > 1) {
                int bendLeft_idx = i * num_width_points + (j - 2);
                springs.emplace_back(&point_masses[idx], &point_masses[bendLeft_idx], BENDING);
            }
            if (i > 1) {
                int bendAbove_idx = (i - 2) * num_width_points + j;
                springs.emplace_back(&point_masses[idx], &point_masses[bendAbove_idx], BENDING);
            }
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.


  // TODO (Part 2): Use Verlet integration to compute new point mass positions


  // TODO (Part 4): Handle self-collisions.


  // TODO (Part 3): Handle collisions with other primitives.


  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  // 2.1: Compute total force acting on each point mass.
  for (auto& pm : point_masses) {
      // Reset force accumulator.
      pm.forces = Vector3D(0, 0, 0);
      // Add external forces F = m * a.
      for (const auto& accel : external_accelerations) {
          pm.forces += mass * accel;
      }
  }

  // Apply spring correction forces.
  for (auto& spring : springs) {
      // Skip spring if the corresponding constraint type is disabled.
      bool enabled = false;
      switch (spring.spring_type) {
      case STRUCTURAL:
          enabled = cp->enable_structural_constraints;
          break;
      case SHEARING:
          enabled = cp->enable_shearing_constraints;
          break;
      case BENDING:
          enabled = cp->enable_bending_constraints;
          break;
      }
      if (!enabled) continue;

      // Compute current spring length and the difference from rest.
      Vector3D vec = spring.pm_a->position - spring.pm_b->position;
      double currLength = vec.norm();
      double ks = cp->ks;
      // Make bending constraints weaker.
      if (spring.spring_type == BENDING)
          ks *= 0.2;
      double diff = currLength - spring.rest_length;

      // Compute the spring force magnitude using Hooke's law.
      double forceMagnitude = ks * diff;
      // Avoid division by zero.
      Vector3D dir = (currLength != 0) ? vec / currLength : Vector3D(0, 0, 0);
      Vector3D force = forceMagnitude * dir;

      // Apply equal and opposite spring forces.
      spring.pm_a->forces -= force;
      spring.pm_b->forces += force;
  }

  // 2.2: Verlet integration to compute new positions.
  for (auto& pm : point_masses) {
      // Do not update positions for pinned masses.
      if (pm.pinned)
          continue;
      Vector3D temp = pm.position; // Store current position to update last_position later.
      Vector3D acceleration = pm.forces / mass;
      // Compute damping factor (cp->damping is a percentage; convert it to a fraction).
      double damping = 1 - cp->damping / 100.0;
      // Verlet integration update: x_new = x_current + (1-damping)*(x_current - x_last) + a * (delta_t)^2.
      pm.position = pm.position + damping * (pm.position - pm.last_position) + acceleration * delta_t * delta_t;
      pm.last_position = temp;
  }

  // 2.3: Constrain the point masses so that spring lengths don't exceed 10% over their rest lengths.
  for (auto& spring : springs) {
      // Skip springs with disabled constraints.
      bool enabled = false;
      switch (spring.spring_type) {
      case STRUCTURAL:
          enabled = cp->enable_structural_constraints;
          break;
      case SHEARING:
          enabled = cp->enable_shearing_constraints;
          break;
      case BENDING:
          enabled = cp->enable_bending_constraints;
          break;
      }
      if (!enabled) continue;

      PointMass* pm_a = spring.pm_a;
      PointMass* pm_b = spring.pm_b;
      Vector3D diff = pm_a->position - pm_b->position;
      double currLength = diff.norm();
      double maxLength = spring.rest_length * 1.1; // Allow 10% stretch.
      if (currLength > maxLength) {
          // Compute the correction distance.
          double correctionDist = currLength - maxLength;
          // Normalized correction direction.
          Vector3D correctionDir = (currLength != 0) ? diff / currLength : Vector3D(0, 0, 0);
          Vector3D correctionVec = correctionDist * correctionDir;

          // Apply correction: if both masses are free, split equally   
          if (!pm_a->pinned && !pm_b->pinned) {
              pm_a->position -= correctionVec * 0.5;
              pm_b->position += correctionVec * 0.5;
          }
          else if (!pm_a->pinned && pm_b->pinned) {
              pm_a->position -= correctionVec;
          }
          else if (pm_a->pinned && !pm_b->pinned) {
              pm_b->position += correctionVec;
          }
          // If both are pinned, no correction is applied.
      }
  }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.

}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.

}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

  return 0.f; 
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}

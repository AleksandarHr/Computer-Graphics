#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 
  Matrix3x3 translation = Matrix3x3::identity();
  translation[2].x = vspan - centerX;
  translation[2].y = vspan - centerY;

  Matrix3x3 scaling = Matrix3x3::identity();
  scaling[0].x = 1 / (2 * vspan);
  scaling[1].y = 1 / (2 * vspan);

  set_svg_2_norm(scaling * translation);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462

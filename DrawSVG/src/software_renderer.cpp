#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;
  transformation_stack.push(svg_2_screen);
  memset(&supersample_target[0], 255, supersample_target.size() * sizeof supersample_target[0]);

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }
  transformation = transformation_stack.top();
  
  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();
  transformation_stack.pop();
}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
    this->sample_rate = sample_rate;
    this->supersample_w = this->target_w * this->sample_rate;
    this->supersample_h = this->target_h * this->sample_rate;
    this->supersample_target.resize(4 * this->supersample_w * this->supersample_h);
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  this->supersample_w = this->target_w * this->sample_rate;
  this->supersample_h = this->target_h * this->sample_rate;
  this->supersample_target.resize(4 * this->supersample_w * this->supersample_h);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
    transformation_stack.push(transformation_stack.top() * element->transform);
    transformation = transformation_stack.top();

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

  transformation_stack.pop();
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  sx *= sample_rate;
  sy *= sample_rate;

  for (int xi = sx; xi < sx + sample_rate; ++xi) {
      for (int yi = sy; yi < sy + sample_rate; ++yi) {
          // fill sample
          size_t pix_position = 4 * (xi + yi * supersample_w);
          supersample_target[pix_position] = (uint8_t)(color.r * 255);
          supersample_target[pix_position + 1] = (uint8_t)(color.g * 255);
          supersample_target[pix_position + 2] = (uint8_t)(color.b * 255);
          supersample_target[pix_position + 3] = (uint8_t)(color.a * 255);
      }
  }

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
    rasterize_line_bresenham(x0, y0, x1, y1, color);
}

void SoftwareRendererImp::rasterize_line_bresenham( float x0, float y0,
                                                    float x1, float y1,
                                                    Color color) {
    
    x0 *= sample_rate;
    x1 *= sample_rate;
    y0 *= sample_rate;
    y1 *= sample_rate;
    
    // We allow lines of slopes 0 <= slope <= 1 only
    bool slopeWithinRange = abs(x1 - x0) >= abs(y1 - y0);
    if (!slopeWithinRange) {
        swap(x0, y0);
        swap(x1, y1);
    }

    // Line should be such that only increases it's x coordinate (e.g. going left-to-right)
    if (x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    int y = y0;
    int epsilon = 0;
    int sign = dy > 0 ? 1 : -1;
    dy = abs(dy);

    for (int x = x0; x <= x1; x++) {
        // Draw point
        if (slopeWithinRange) {
            fill_sample(x, y, color);
        }
        else {
            fill_sample(y, x, color);
        }

        // Update coordinates for next point
        epsilon += dy;
        if ((epsilon << 1) >= dx) {
            y += sign;
            epsilon -= dx;
        }
    }
}

void SoftwareRendererImp::fill_sample(int x, int y, Color color) {

    // check bounds
    if (x < 0 || x >= this->supersample_w) return;
    if (y < 0 || y >= this->supersample_h) return;

    size_t pos = 4 * (x + y * this->supersample_w);
    Color from = color;

    Color src, to;
    src.r = supersample_target[pos] / 255.0f;
    src.g = supersample_target[pos + 1] / 255.0f;
    src.b = supersample_target[pos + 2] / 255.0f;
    src.a = supersample_target[pos + 3] / 255.0f;

    to.r = (1.0f - from.a) * src.r + from.r * from.a;
    to.g = (1.0f - from.a) * src.g + from.g * from.a;
    to.b = (1.0f - from.a) * src.b + from.b * from.a;
    to.a = 1.0f - (1.0f - from.a) * (1.0f - src.a);

    supersample_target[pos] = (uint8_t)(to.r * 255);
    supersample_target[pos + 1] = (uint8_t)(to.g * 255);
    supersample_target[pos + 2] = (uint8_t)(to.b * 255);
    supersample_target[pos + 3] = (uint8_t)(to.a * 255);
}

// TODO: Improve quality
void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    x0 *= sample_rate;
    y0 *= sample_rate;
    x1 *= sample_rate;
    y1 *= sample_rate;
    x2 *= sample_rate;
    y2 *= sample_rate;
    
    // Find min/max coordinates to determine bounding box of the triangle
    float maxX = ceil(max(max(x0, x1), x2));
    float maxY = ceil(max(max(y0, y1), y2));
    float minX = floor(min(min(x0, x1), x2)) + 0.5f;
    float minY = floor(min(min(y0, y1), y2)) + 0.5f;

    // Iterate over each sample and determine if it is within the triangle specified
    // TODO: Improve algorithm -- 'early out' and 'early in' boxes
    for (float y = minY; y <= maxY; y+=1) {
        for (float x = minX; x <= maxX; x+=1) {
            int turns = point_plane_check(x0, y0, x1, y1, x, y) +
                point_plane_check(x1, y1, x2, y2, x, y) +
                point_plane_check(x2, y2, x0, y0, x, y);

            // If it's within the triangle, color
            bool within_triangle = (abs(turns) == 3);
            if (within_triangle) {
                fill_sample(x, y, color);
            }
        }
    }
}

// Check which side the point (x2,y2) is from the line (x0,y0) - (x1,y1)
int SoftwareRendererImp::point_plane_check(float x0, float y0,
                                             float x1, float y1,
                                             float x2, float y2) {
    float dx01 = x1 - x0;
    float dy01 = y1 - y0;
    float dx02 = x2 - x0;
    float dy02 = y2 - y0;
    float cross_product = (dx01 * dy02) - (dy01 * dx02);
    return cross_product > 0 ? 1 : -1;
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling

    size_t num_samples = sample_rate * sample_rate;
    for (size_t x = 0; x <= supersample_w - sample_rate; x += sample_rate) {
        for (size_t y = 0; y <= supersample_h - sample_rate; y += sample_rate) {
            
            // Iterate over the sample_rate x sample_rate box
            uint16_t r = 0, g = 0, b = 0, a = 0;
            for (size_t xi = 0; xi < sample_rate; ++xi) {
                for (size_t yi = 0; yi < sample_rate; ++yi) {
                    size_t sample_position = 4 * (x + xi + (y + yi) * supersample_w);
                    r += supersample_target[sample_position];
                    g += supersample_target[sample_position + 1];
                    b += supersample_target[sample_position + 2];
                    a += supersample_target[sample_position + 3];
                }
            }

            // Compute average of r,g,b,a channels over the sampling box
            r = r / num_samples;
            g = g / num_samples;
            b = b / num_samples;
            a = a / num_samples;

            // Map the computed value to the relevant original pixel
            size_t pix_x = x / sample_rate;
            size_t pix_y = y / sample_rate;
            size_t pix_position = 4 * (pix_x + pix_y * target_w);
            render_target[pix_position] = (uint8_t)(r);
            render_target[pix_position + 1] = (uint8_t)(g);
            render_target[pix_position + 2] = (uint8_t)(b);
            render_target[pix_position + 3] = (uint8_t)(a);
        }
    }
}


} // namespace CMU462

#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // Task 7: Implement this
  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // Fill in mipmaps with appropriate data
  // MimMap[0] is the original image
  for (size_t idx = 1; idx < tex.mipmap.size(); idx++) {
      MipLevel& prev = tex.mipmap[idx - 1];
      MipLevel& curr = tex.mipmap[idx];

      for (size_t x = 0; x < curr.width; x++) {
          // corresponding x-coordinate in previous mipmap
          size_t prev_x = 2 * x;
          for (size_t y = 0; y < curr.height; y++) {
              // corresponding y-coordinate in previous mipmap
              size_t prev_y = 2 * y;
              uint16_t r = 0, g = 0, b = 0, a = 0;
              // Sum the values of the four corresponding pixels in a 2x2 window
              for (int i = 0; i < 2; i++) {
                  for (int j = 0; j < 2; j++) {
                      size_t prev_pix_pos = 4 * (prev_x + i + (prev_y + j) * prev.width);
                      r += prev.texels[prev_pix_pos];
                      g += prev.texels[prev_pix_pos + 1];
                      b += prev.texels[prev_pix_pos + 2];
                      a += prev.texels[prev_pix_pos + 3];
                  }
              }
              // Update the current mipmap pixel values with averaged values
              size_t curr_pix_pos = 4 * (x + y * curr.width);
              curr.texels[curr_pix_pos] = r / 4;
              curr.texels[curr_pix_pos + 1] = g / 4;
              curr.texels[curr_pix_pos + 2] = b / 4;
              curr.texels[curr_pix_pos + 3] = a / 4;
          }
      }
  }
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
    if (level > tex.mipmap.size()) {
        // return magenta for invalid level
        return Color(1, 0, 1, 1);
    }

    MipLevel& mip = tex.mipmap[level];
    float adjusted_u = u * mip.width;
    float adjusted_v = v * mip.height;
    size_t texel_x = adjusted_u < 0.5f ? 0 : round(adjusted_u - 0.5f);
    size_t texel_y = adjusted_v < 0.5f ? 0 : round(adjusted_v - 0.5f);
    Color c;
    uint8_to_float(&c.r, &mip.texels[4 * (texel_x + mip.width * texel_y)]);
    return c;
}

Color Sampler2DImp::sample_trilinear(Texture& tex,
    float u, float v,
    float u_scale, float v_scale) {

    // compute mipmap level given u_scale and v_scale parameters values
    float d = max(log2f(max(u_scale, v_scale)), 0.0f);

//    float d = max(0.f, log2f(max(u_scale, v_scale)));
    // two mipmap levels to be interpolated between
    int first_level = floor(d);
    int second_level = first_level + 1;

    // interpolation using the fractional part of the mipmap level
    float w = d - first_level;

    Color c0 = sample_bilinear(tex, u, v, first_level);
    Color c1 = sample_bilinear(tex, u, v, second_level);

    return (1 - w) * c0 + w * c1;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering

  // return magenta for invalid level
    if (level > tex.mipmap.size()) {
        // return magenta for invalid level
        return Color(1, 0, 1, 1);
    }

    MipLevel& mip = tex.mipmap[level];
    float texel_u = u * mip.width;
    float texel_v = v * mip.height;
    size_t texel_i = texel_u < 0.5f ? 0 : floor(texel_u - 0.5f);
    size_t texel_j = texel_v < 0.5f ? 0 : floor(texel_v - 0.5f);
    float texel_s = texel_u - texel_i - 0.5f;
    float texel_t = texel_v - texel_j - 0.5f;

    Color c00, c01, c10, c11;
    uint8_to_float(&c00.r, &mip.texels[4 * (texel_i + mip.width * texel_j)]);
    uint8_to_float(&c01.r, &mip.texels[4 * (texel_i + 1 + mip.width * texel_j)]);
    uint8_to_float(&c10.r, &mip.texels[4 * (texel_i + mip.width * (texel_j + 1))]);
    uint8_to_float(&c11.r, &mip.texels[4 * (texel_i + 1 + mip.width * (texel_j + 1))]);

    Color f = (1 - texel_t)*((1 - texel_s) * c00 + texel_s * c01) + 
        texel_t * ((1-texel_s)*c10 + texel_s * c11);
    return f;
}

} // namespace CMU462

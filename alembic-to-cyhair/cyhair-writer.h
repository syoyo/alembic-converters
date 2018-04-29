/*
MIT License

Copyright (c) 2018 Syoyo Fujita

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef CYHAIR_WRITER_H_
#define CYHAIR_WRITER_H_

#include <cstdint>  // C++11
#include <cstdlib>
#include <vector>

namespace cyhair_writer {

// Cyhair
#define CYHAIR_INFO_HEADER_SIZE (88)

struct CyhairHeader {
  char signagure[4];  // "HAIR"
  uint32_t hair_count;
  uint32_t point_count;
  uint32_t flags;

  // default values
  uint32_t d_segments;
  float d_thickness;
  float d_transparency;
  float d_color[3];

  char info[CYHAIR_INFO_HEADER_SIZE];  // infos. pad to 128bytes.
};

///
/// Save curves as CyHair format.
/// TODO(syoyo): Write per-CV thickness, transparency and color.
///
bool SaveAsCyhair(const std::string &filename, std::vector<float> points,
                  std::vector<uint32_t> num_points_per_segments);

}  // namespace cyhair_writer

#endif  // CYHAIR_WRITER_H_

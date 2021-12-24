#include <math.h>
#include <ostream>
#include <vector>
#include <iostream>

#include "Magick++.h"
#include "Magick++/Functions.h"
#include "MagickCore/image.h"

#include "sample_entropy_calculator2d.h"
#include "utils.h"

struct Argument {
  unsigned m;
  double r;
  unsigned x;
  unsigned y;
  unsigned width;
  unsigned height;
  unsigned moving_step_size;
  unsigned dilation_factor;
  std::string image_filename;
  sampen::OutputLevel output_level;
} arg;


char usage[] = \
" --input <INPUT_IMAGE_FILENAME> -r <R> (default = 0.3)\\\n"
"    -m <M> (default = 1) --moving-step-size <STEP_SIZE>(default = 1)\\\n"
"    --dilation-factor <FACTOR> (default = 1)\n\n"
"Arguments:\n"
"    --output-level <LEVEL> (default = 0)\n"
"        Can be 0, 1 or 2, amongst which 0 is most silent and 2 is most verbose.\n"
"    --input <IMAGE_FILENAME>\n"
"        The filename of the image to calculate SampEn2D.\n"
"    -r <R> (default = 0.3)\n"
"        Threshold for template matching. Note that in actual computation,\n"
"        this value will be scaled by the standard deviation of the input image.\n"
"    -m <M> (default = 1)\n"
"        Template length. The size sliding window will be mxm and (m+1)x(m+1).\n"
"    --moving-step-size <STEP_SIZE> (default = 1)\n"
"        The number of pixels that the sliding windows move at one time.\n"
"    -x <X> -y <Y> -w <W> -h <H>\n"
"        The position and size of the rectangle to crop the original image, and\n"
"        the SampEn2D will be calculated on the cropped image. If these values\n"
"        are not presented, the whole image will be used.\n"
"    --dilation-factor <FACTOR> (default = 1)\n"
"        This implementation allows holes in the sliding window, as like dilated\n"
"        convolution in deep learning.\n\n"
"Options:\n"
"    --help\n"
"        Display this message.\n";

void ParseArgument(int argc, char *argv[]) {
  sampen::ArgumentParser parser(argc, argv);
  if (parser.isOption("--help")) {
    std::cout << "Usage: " << argv[0] << usage;
    exit(0);
  }
  arg.m = parser.getArgLong("-m", 1);
  arg.r = parser.getArgDouble("-r", 0.3);
  arg.x = parser.getArgLong("-x", 0);
  arg.y = parser.getArgLong("-y", 0);
  arg.width = parser.getArgLong("-w", 0);
  arg.height = parser.getArgLong("-h", 0);
  arg.moving_step_size = parser.getArgLong("--moving-step-size", 1);
  arg.dilation_factor = parser.getArgLong("--dilation-factor", 1);
  arg.image_filename = parser.getArg("--input");
  arg.output_level = static_cast<sampen::OutputLevel>(
      parser.getArgLong("--output-level", 0));
  if (arg.image_filename.empty()) {
    std::cerr << "Input image file name is required.\n";
    std::cerr << "Usage: " << argv[0] << usage;
    exit(-1);
  }
}

int main(int argc, char *argv[]) {
  Magick::InitializeMagick(argv[0]);
  Magick::Image image;
  ParseArgument(argc, argv);

  image.read(arg.image_filename);

  unsigned width = image.columns();
  unsigned height = image.rows();
  
  image.type(Magick::GrayscaleType);

  if (arg.width == 0) arg.width = width;
  if (arg.height == 0) arg.height = height;
  if (arg.x + arg.width > width || arg.y + arg.height > height) {
    MSG_ERROR(-1, "Specified rectangle is outside the image geometry.\n");
  } 
  // Print arguments.
  std::cout << "========================================"
      << "========================================\n";
  std::cout << "Sample Entropy (2D) Computation Setting:\n";
  std::cout << "\tfilename: " << image.fileName() << std::endl;
  std::cout << "\twidth: " << image.columns() << std::endl;
  std::cout << "\theight: " << image.rows() << std::endl;
  std::cout << "\tr: " << arg.r << std::endl;
  std::cout << "\tm: " << arg.m << std::endl;
  std::cout << "\tx: " << arg.x << std::endl;
  std::cout << "\ty: " << arg.y << std::endl;
  std::cout << "\tw: " << arg.width << std::endl;
  std::cout << "\th: " << arg.height << std::endl;
  std::cout << "\tmoving-step-size: " << arg.moving_step_size << std::endl;
  std::cout << "\tdilation-factor: " << arg.dilation_factor << std::endl;

  auto num_channels = image.channels();
  auto pixels = image.getConstPixels(arg.x, arg.y, arg.width, arg.height);
  std::vector<int> data(arg.width * arg.height * num_channels);

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = static_cast<int>(pixels[i]);
  }
  auto var = sampen::ComputeVariance(data);
  int r = sqrt(var) * arg.r;
  sampen::SampleEntropyCalculator2DDirect<int> sec2dd(
      data.begin(), data.end(), arg.m, r, width, height, arg.moving_step_size,
      arg.dilation_factor, arg.output_level);
  sec2dd.ComputeSampleEntropy();
  std::cout << sec2dd.get_result_str();
  std::cout << "========================================"
      << "========================================\n";
}
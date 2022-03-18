#include <iostream>
#include <math.h>
#include <memory>
#include <ostream>
#include <fstream>
#include <vector>

#include "Magick++.h"
#include "Magick++/Functions.h"
#include "Magick++/Geometry.h"
#include "Magick++/Image.h"
#include "MagickCore/image.h"

#include "MagickCore/pixel.h"
#include "sample_entropy_calculator2d.h"
#include "utils.h"

using namespace sampen;

struct Argument {
  unsigned m;
  double r;
  unsigned x;
  unsigned y;
  unsigned width;
  unsigned height;
  unsigned moving_step_size;
  unsigned dilation_factor;
  bool sampling1;
  bool sampling2;
  unsigned multiscale_depth;
  double multiscale_factor;
  unsigned sample_size;
  unsigned sample_num;
  std::string image_filename;
  sampen::OutputLevel output_level;
} arg;

// clang-format off
char usage[] =
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
    "    --sample-size <N0>\n"
    "        The number of templates drawn from template set to calculate SampEn2D on.\n"
    "    --sample-num <N1>\n"
    "        The number of experiments conducted for sampling based methods.\n"
    "    --dilation-factor <FACTOR> (default = 1)\n"
    "        This implementation allows holes in the sliding window, as like dilated\n"
    "        convolution in deep learning.\n\n"
    "Options:\n"
    "    --sampling1\n"
    "        Use sampling method (approach 1).\n"
    "    --sampling2\n"
    "        Use sampling method (approach 2).\n"
    "    --help\n"
    "        Display this message.\n";
// clang-format on

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
  arg.output_level =
      static_cast<sampen::OutputLevel>(parser.getArgLong("--output-level", 0));
  arg.sampling1 = parser.isOption("--sampling1");
  arg.sampling2 = parser.isOption("--sampling2");
  arg.multiscale_depth = parser.getArgLong("--multiscale-depth", 1);
  arg.multiscale_factor = parser.getArgDouble("--multiscale-factor", 0.5);
  if (arg.multiscale_depth <= 0) {
    MSG_ERROR(-1, "Invalid argument for --multiscale-depth, positive integer "
                  "required.\n");
  }
  if (arg.image_filename.empty()) {
    std::cerr << "Input image file name is required.\n";
    std::cerr << "Usage: " << argv[0] << usage;
    exit(-1);
  }
  if (arg.sampling1 || arg.sampling2) {
    arg.sample_size = parser.getArgLong("--sample-size", 1024);
    arg.sample_num = parser.getArgLong("--sample-num", 20);
  }
}

std::vector<double> ComputeMultiscaleSampEn2D(Magick::Image image, double r,
                                              int m, int depth, double factor,
                                              bool sampling) {
  image.type(Magick::GrayscaleType);
  int width = image.columns();
  int height = image.rows();
  std::vector<double> result;
  PrintSeperator('-');
  for (int i = 0; i < depth; ++i) {
    image.resize(Magick::Geometry(width, height));
    width = image.columns();
    height = image.rows();
    if (height <= m + 2 || width <= m + 2) {
        break;
    }
    auto pixels = image.getConstPixels(0, 0, width, height);
    int num_channels = image.channels();
    std::vector<int> data(width * height * 1);
    for (size_t i = 0; i < data.size(); ++i) {
      data[i] = static_cast<int>(pixels[i * num_channels]);
    }

    std::shared_ptr<sampen::SampleEntropyCalculator2D<int>> calculator;
    if (sampling) {
      calculator =
          std::make_shared<SampleEntropyCalculator2DSamplingDirect<int>>(
              data.cbegin(), data.cend(), m, r, width, height, 1, 1,
              arg.sample_size, arg.sample_num, 0., 0., 0., arg.output_level);
    } else {
      calculator = std::make_shared<SampleEntropyCalculator2DDirect<int>>(
          data.cbegin(), data.cend(), m, r, width, height, 1, 1,
          arg.output_level);
    }
    std::cout << "SampEn2D (Scale " << i
        << ", h: " << height << ", w: " << width << "): "
        << calculator->get_entropy()
        << std::endl;
    result.push_back(calculator->get_entropy());
    if (factor - 0.5 < 1e-8 && factor - 0.5 > -1e-8) {
        width /= 2;
        height /= 2;
    } else {
        width = static_cast<int>(static_cast<double>(width) * factor);
        height = static_cast<int>(static_cast<double>(height) * factor);
    }
  }
  return result;
}

int main(int argc, char *argv[]) {
  Magick::InitializeMagick(argv[0]);
  Magick::Image image;
  ParseArgument(argc, argv);

  image.read(arg.image_filename);

  unsigned width = image.columns();
  unsigned height = image.rows();

  image.type(Magick::GrayscaleType);
  image.channel(MagickCore::RedChannel);

  if (arg.width == 0)
    arg.width = width;
  if (arg.height == 0)
    arg.height = height;
  if (arg.x + arg.width > width || arg.y + arg.height > height) {
    MSG_ERROR(-1, "Specified rectangle is outside the image geometry.\n");
  }
  // Print arguments.
  PrintSeperator('=');
  std::cout << "Sample Entropy (2D) Computation Setting:\n";
  std::cout << "\tfilename: " << image.fileName() << std::endl;
  std::cout << "\twidth: " << image.columns() << std::endl;
  std::cout << "\theight: " << image.rows() << std::endl;
  std::cout << "\tnum_channels: " << image.channels() << std::endl;
  std::cout << "\tr: " << arg.r << std::endl;
  std::cout << "\tm: " << arg.m << std::endl;
  std::cout << "\tx: " << arg.x << std::endl;
  std::cout << "\ty: " << arg.y << std::endl;
  std::cout << "\tw: " << arg.width << std::endl;
  std::cout << "\th: " << arg.height << std::endl;
  std::cout << "\tmoving-step-size: " << arg.moving_step_size << std::endl;
  std::cout << "\tdilation-factor: " << arg.dilation_factor << std::endl;
  std::cout << "\tmultiscale-depth: " << arg.multiscale_depth << std::endl;
  if (arg.multiscale_depth) {
      std::cout << "\tmultiscale-factor: " << arg.multiscale_factor << std::endl;
  }


  auto num_channels = image.channels();
  auto pixels = image.getConstPixels(arg.x, arg.y, arg.width, arg.height);
  std::vector<int> data(arg.width * arg.height * 1);
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = static_cast<int>(pixels[i * num_channels]);
  }
  auto var = sampen::ComputeVariance(data);
  std::cout << "\n\tstandard deviation: " << sqrt(var) << std::endl;
  arg.r = sqrt(var) * arg.r;


  if (arg.sampling1) {
    sampen::SampleEntropyCalculator2DDirect<int> sec2dd(
        data.begin(), data.end(), arg.m, arg.r, arg.width, arg.height,
        arg.moving_step_size, arg.dilation_factor, arg.output_level);
    sec2dd.ComputeSampleEntropy();
    std::cout << sec2dd.get_result_str();

    double sampen2d = sec2dd.get_entropy();
    double a_norm = sec2dd.get_a_norm();
    double b_norm = sec2dd.get_b_norm();
    sampen::SampleEntropyCalculator2DSamplingDirect<int> sec2dds(
        data.cbegin(), data.cend(), arg.m, arg.r, arg.width, arg.height,
        arg.moving_step_size, arg.dilation_factor, arg.sample_size,
        arg.sample_num, sampen2d, a_norm, b_norm, arg.output_level);
    std::cout << sec2dds.get_result_str();
  }
  ComputeMultiscaleSampEn2D(image, arg.r, arg.m, arg.multiscale_depth,
                            arg.multiscale_factor, false);
  PrintSeperator('=');
}

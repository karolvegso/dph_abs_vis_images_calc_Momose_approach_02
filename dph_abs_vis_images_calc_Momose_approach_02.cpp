// dph_abs_vis_images_calc_Momose_approach_02.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <chrono>

using namespace std;

int main()
{
    // define path to folder with all foreground data or all subfolders
    string path_to_fg_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/pp/fg/");
    // define path to folder with all background data or all subfolders
    string path_to_bg_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/bg/");
    // print path to folder with all foreground folders
    std::cout << path_to_fg_folder << "\n";
    // print path to folder with all background folders
    std::cout << path_to_bg_folder << "\n";

    // define path to output folder with output differential phase (dph) images
    string path_to_output_folder("d:/XTI_Momose_lab/BL28B2_2017A/sort_data/pp_dph_abs_vis_Momose/");

    // output image name root - differential phase image or dph image
    string image_output_dph_name_root = "dph";
    // output image name root - absorption image or abs image
    string image_output_abs_name_root = "abs";
    // output image name root - visibility image or vis image
    string image_output_vis_name_root = "vis";
    // final ouput image name - differential phase image or dph image
    string image_output_dph_name;
    // final ouput image name - absorption image or abs image
    string image_output_abs_name;
    // final ouput image name - visibility image or vis image
    string image_output_vis_name;
    // extension of output image
    string image_output_extension = ".raw";

    // define size of the raw unsigned 16 bit images
    const unsigned int no_cols = 1536; // in pixels, in horizontal direction
    const unsigned int no_rows = 512; // in pixels, in vertical direction
    // total number of pixels in single image
    const unsigned int no_pixels = no_cols * no_rows;
    // total number of bytes in single image, we consider 16 bit values per pixel = 2 bytes
    const unsigned int no_bytes = 2 * no_pixels;

    // define number of initial and final subfolder for foreground
    unsigned int no_subfolder_fg_initial = 1;
    unsigned int no_subfolder_fg_final = 11200;

    // define number of initial and final folder for background
    unsigned int no_subfolder_bg_initial = 1;
    unsigned int no_subfolder_bg_final = 11200;

    // number of digits in subfolder name for foreground
    string::size_type no_subfolder_digits_fg = 6;
    // number of digits in subfolder name for background
    string::size_type no_subfolder_digits_bg = 6;

    // number of steps in fringe scanning technique
    const unsigned int M = 5;

    // calculate differential phase image for foreground
    // fringe scanning defined from initial value
    const unsigned int M_fg_initial = 1;
    // fringe scanning defined to final value
    const unsigned int M_fg_final = M;
    // number of steps in fringe scanning for foreground
    const unsigned int M_fg = M_fg_final - M_fg_initial + 1;

    // calculate differential phase image for background
    // fringe scanning defined from initial value
    const unsigned int M_bg_initial = 1;
    // fringe scanning defined to final value
    const unsigned int M_bg_final = M;
    // number of steps in fringe scanning for background
    const unsigned int M_bg = M_bg_final - M_bg_initial + 1;

    // define root name of images for foreground
    string root_image_name_fg("a");
    // define root name of images for background
    string root_image_name_bg("a");

    // number of digits in image name for foreground
    string::size_type no_image_digits_fg = 6;
    // number of digits in image name for background
    string::size_type no_image_digits_bg = 6;

    // define image extensions
    // image extension for foreground
    string image_extension_fg = ".raw";
    // image extension for background
    string image_extension_bg = ".raw";

    // allocate image buffer for foreground
    auto image_buffer_fg = new unsigned short int[no_pixels][M_fg];
    // allocate image buffer for background
    auto image_buffer_bg = new unsigned short int[no_pixels][M_fg];

    // allocate phase buffer for foreground
    double* phase_buffer_fg = new double[no_pixels];
    // allocate phase buffer for background
    double* phase_buffer_bg = new double[no_pixels];
    // allocate amplitude buffer for foreground
    double* amp_buffer_fg = new double[no_pixels];
    // allocate amplitude buffer for background
    double* amp_buffer_bg = new double[no_pixels];
    // allocate offset buffer for foreground
    double* offset_buffer_fg = new double[no_pixels];
    // allocate offset buffer for background
    double* offset_buffer_bg = new double[no_pixels];

    // allocate memory for differential phase image
    double* dph_image = new double[no_pixels];
    // allocate memory for absorption image
    double* abs_image = new double[no_pixels];
    // allocate memory for visibility image
    double* vis_image = new double[no_pixels];

    // define phase for foreground
    double phase_step_fg = (2 * M_PI) / M_fg;
    // define phase_step for background
    double phase_step_bg = (2 * M_PI) / M_bg;

    // auxiliary variables for iteration through subfolder name for foreground
    string subfolder_name(no_subfolder_digits_fg, '0');
    string subfolder_number = "";
    string::size_type counter_digits = 0;
    string::size_type difference = 0;
    string::size_type counter = 0;
    string path_to_fg_subfolder = "";

    // auxiliary variables for iteration through M_fg images
    int counter_image = 0;
    string image_name = root_image_name_fg;
    string image_name_number(no_image_digits_fg, '0');
    string image_number = "";
    // counter_digits, difference and counter variables are taken from iterations through subfolders
    string path_to_fg_image = "";

    // auxiliary variables for iteration through subfolder name for background
    string path_to_bg_subfolder = "";

    // auxiliary variables for iteration through M_bg images
    image_name = root_image_name_bg;
    image_name_number = string(no_image_digits_bg, '0');
    image_number = "";
    // counter_digits, difference and counter variables are taken from iterations through subfolders
    string path_to_bg_image = "";

    // declare auxiliary variable for output image
    string path_to_output_image = "";

    // declare auxiliary variable for output dph image
    string path_to_output_dph_image = "";
    // declare auxiliary variable for output abs image
    string path_to_output_abs_image = "";
    // declare auxiliary variable for output vis image
    string path_to_output_vis_image = "";

    // go through all foreground subfolders and foreground images
    for (unsigned int index_0 = no_subfolder_fg_initial; index_0 <= no_subfolder_fg_final; index_0++) {
        // start to measure elapsed time at the beginning
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // initialize subfolder name to "000000"
        subfolder_name = string(no_subfolder_digits_fg, '0');
        // typcast integer value to string, convert integer value to string
        subfolder_number = std::to_string(index_0);
        // initialize digits counter
        counter_digits = subfolder_number.size();
        // initialize difference
        difference = no_subfolder_digits_fg - counter_digits;
        // initialize counter
        counter = 0;
        // generate subfolder name
        for (string::size_type index_1 = difference; index_1 < no_subfolder_digits_fg; index_1++) {
            subfolder_name[index_1] = subfolder_number[counter];
            counter++;
        }
        // generate path to foreground subfolder 
        path_to_fg_subfolder = path_to_fg_folder + subfolder_name;
        // print final path to foreground subfolder
        std::cout << path_to_fg_subfolder << "\n";
        // initilize counter for one of M images
        counter_image = 0;
        // open images in foreground subfolder
        for (unsigned int index_2 = M_fg_initial; index_2 <= M_fg_final; index_2++) {
            // initialize image name to root value
            image_name = root_image_name_fg;
            // initialize image number to "000000"
            image_name_number = string(no_image_digits_fg, '0');
            // typcast integer value to string, convert integer value to string
            image_number = std::to_string(index_2);
            // initialize digits counter
            counter_digits = image_number.size();
            // initialize difference
            difference = no_image_digits_fg - counter_digits;
            // initialize counter
            counter = 0;
            // generate image name number
            for (string::size_type index_3 = difference; index_3 < no_image_digits_fg; index_3++) {
                image_name_number[index_3] = image_number[counter];
                counter++;
            }
            // concatenate root image value and image number
            image_name += image_name_number;
            // generate path to foreground image 
            path_to_fg_image = path_to_fg_subfolder + "/" + image_name + image_extension_fg;
            // print path to image for foregrouund
            //std::cout << path_to_fg_image << "\n";
            //************************************************************************************
            // read binary images
            //************************************************************************************
            ifstream raw_image(path_to_fg_image, ios::out | ios::binary);
            /*streampos begin, end;
            begin = raw_image.tellg();
            raw_image.seekg(0, ios::end);
            end = raw_image.tellg();*/
            //std::cout << "Size of the raw image is: " << (end - begin) << " bytes.\n";
            if (raw_image.is_open())
            {
                //unsigned int counter_pixel = 0;
                //raw_image.seekg(0, ios::beg);
                //while (raw_image.read(reinterpret_cast<char*>(&image_buffer_fg[counter_pixel][counter_image]), sizeof(uint16_t))) { // Read 16-bit integer values from file
                //    counter_pixel++;
                //}
                //raw_image.close();
                for (unsigned int counter_pixel = 0; counter_pixel < no_pixels; counter_pixel++) {
                    raw_image.read((char*)&image_buffer_fg[counter_pixel][counter_image], sizeof(unsigned short int));
                }
                raw_image.close();
            }
            else {
                std::cout << "Warning: Unable to open raw image file!!!" << "\n";
            }
            //************************************************************************************
            // end of reading of binary images
            //************************************************************************************
            // increase image counter by one
            counter_image++;
        }

        // go through background subfolder and background images
        // initialize subfolder name to "000000"
        subfolder_name = string(no_subfolder_digits_bg, '0');
        // typcast integer value to string, convert integer value to string
        subfolder_number = std::to_string(index_0);
        // initialize digits counter
        counter_digits = subfolder_number.size();
        // initialize difference
        difference = no_subfolder_digits_bg - counter_digits;
        // initialize counter
        counter = 0;
        // generate subfolder name
        for (string::size_type index_1 = difference; index_1 < no_subfolder_digits_bg; index_1++) {
            subfolder_name[index_1] = subfolder_number[counter];
            counter++;
        }
        // generate path to background subfolder 
        path_to_bg_subfolder = path_to_bg_folder + subfolder_name;
        // print final path to background subfolder
        std::cout << path_to_bg_subfolder << "\n";
        // initilize counter for one of M images
        counter_image = 0;
        // open images in background subfolder
        for (unsigned int index_2 = M_bg_initial; index_2 <= M_bg_final; index_2++) {
            // initialize image name to root value
            image_name = root_image_name_bg;
            // initialize image number to "000000"
            image_name_number = string(no_image_digits_bg, '0');
            // typcast integer value to string, convert integer value to string
            image_number = std::to_string(index_2);
            // initialize digits counter
            counter_digits = image_number.size();
            // initialize difference
            difference = no_image_digits_bg - counter_digits;
            // initialize counter
            counter = 0;
            // generate image name number
            for (string::size_type index_3 = difference; index_3 < no_image_digits_bg; index_3++) {
                image_name_number[index_3] = image_number[counter];
                counter++;
            }
            // concatenate root image value and image number
            image_name += image_name_number;
            // generate path to background image 
            path_to_bg_image = path_to_bg_subfolder + "/" + image_name + image_extension_bg;
            // print path to image for background
            //std::cout << path_to_bg_image << "\n";
            //************************************************************************************
            // read binary images
            //************************************************************************************
            ifstream raw_image(path_to_bg_image, ios::out | ios::binary);
            /*streampos begin, end;
            begin = raw_image.tellg();
            raw_image.seekg(0, ios::end);
            end = raw_image.tellg();*/
            //std::cout << "Size of the raw image is: " << (end - begin) << " bytes.\n";
            if (raw_image.is_open())
            {
                //unsigned int counter_pixel = 0;
                //raw_image.seekg(0, ios::beg);
                //while (raw_image.read(reinterpret_cast<char*>(&image_buffer_bg[counter_pixel][counter_image]), sizeof(uint16_t))) { // Read 16-bit integer values from file
                //    counter_pixel++;
                //}
                //raw_image.close();
                raw_image.seekg(0, ios::beg);
                for (unsigned int counter_pixel = 0; counter_pixel < no_pixels; counter_pixel++) {
                    raw_image.read((char*)&image_buffer_bg[counter_pixel][counter_image], sizeof(unsigned short int));
                }
                raw_image.close();
            }
            else {
                std::cout << "Warning: Unable to open raw image file!!!" << "\n";
            }
            //************************************************************************************
            // end of reading of binary images
            //************************************************************************************
            // increase image counter by one
            counter_image++;
        }

        // calculate phase image for foreground
        // initialize complex number z
        std::complex<double> z = 0.0 + 0.0i;
        // initialize complex number z1
        std::complex<double> z1 = 0.0 + 0.0i;
        double iterator = 0.0f;
        double intensity_value = 0.0f;
        double sum_intensity = 0.0f;
        for (unsigned int index_8 = 0; index_8 < no_pixels; index_8++) {
            z = 0.0 + 0.0i;
            sum_intensity = 0.0f;
            for (unsigned int index_9 = 0; index_9 < M_fg; index_9++) {
                iterator = double(index_9);
                intensity_value = double(image_buffer_fg[index_8][index_9]);
                sum_intensity += intensity_value;
                z1 = intensity_value * std::exp(1i * iterator * phase_step_fg);;
                z = z + z1;
            }
            phase_buffer_fg[index_8] = std::arg(z);
            amp_buffer_fg[index_8] = 2.0f * std::abs(z) / double(M_fg);
            offset_buffer_fg[index_8] = sum_intensity / double(M_fg);
        }
        // calculate phase image for background
        // initialize complex number z
        z = 0.0 + 0.0i;
        // initialize complex number z1
        z1 = 0.0 + 0.0i;
        iterator = 0.0f;
        intensity_value = 0.0f;
        sum_intensity = 0.0f;
        for (unsigned int index_8 = 0; index_8 < no_pixels; index_8++) {
            z = 0.0 + 0.0i;
            sum_intensity = 0.0f;
            for (unsigned int index_9 = 0; index_9 < M_bg; index_9++) {
                iterator = double(index_9);
                intensity_value = double(image_buffer_bg[index_8][index_9]);
                sum_intensity += intensity_value;
                z1 = intensity_value * std::exp(1i * iterator * phase_step_bg);
                z = z + z1;
            }
            phase_buffer_bg[index_8] = std::arg(z);
            amp_buffer_bg[index_8] = 2.0f * std::abs(z) / double(M_bg);
            offset_buffer_bg[index_8] = sum_intensity / double(M_bg);
        }

        // calculate differential phase image or dph image
        for (unsigned int index_10 = 0; index_10 < no_pixels; index_10++) {
            dph_image[index_10] = phase_buffer_fg[index_10] - phase_buffer_bg[index_10];
            abs_image[index_10] = offset_buffer_fg[index_10] / offset_buffer_bg[index_10];
            vis_image[index_10] = (amp_buffer_fg[index_10] / offset_buffer_fg[index_10]) / (amp_buffer_bg[index_10] / offset_buffer_bg[index_10]);
        }

        // define name for output dph image for current subfolder
        image_output_dph_name = image_output_dph_name_root + "_" + subfolder_name + image_output_extension;
        // define name for output abs image for current subfolder
        image_output_abs_name = image_output_abs_name_root + "_" + subfolder_name + image_output_extension;
        // define name for output vis image for current subfolder
        image_output_vis_name = image_output_vis_name_root + "_" + subfolder_name + image_output_extension;
        // define path to the output dph image
        path_to_output_dph_image = path_to_output_folder + image_output_dph_name;
        // define path to the output abs image
        path_to_output_abs_image = path_to_output_folder + image_output_abs_name;
        // define path to the output vis image
        path_to_output_vis_image = path_to_output_folder + image_output_vis_name;
        // write differential phase (dph) image
        // set for output, binary data, trunc
        fstream output_dph_image(path_to_output_dph_image, ios::out | ios::binary | ios::trunc);
        if (output_dph_image.is_open())
        {
            // set pointer to the beginning of the image
            output_dph_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_dph_image.write((char*)&dph_image[index_11], sizeof(double));
            }
            output_dph_image.close();
        }
        else {
            std::cout << "Warning: Unable to open dph image file!!!" << "\n";
        }
        // write absorption (abs) image
        // set for output, binary data, trunc
        fstream output_abs_image(path_to_output_abs_image, ios::out | ios::binary | ios::trunc);
        if (output_abs_image.is_open())
        {
            // set pointer to the beginning of the image
            output_abs_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_abs_image.write((char*)&abs_image[index_11], sizeof(double));
            }
            output_abs_image.close();
        }
        else {
            std::cout << "Warning: Unable to open abs image file!!!" << "\n";
        }
        // write visibility (vis) image
        // set for output, binary data, trunc
        fstream output_vis_image(path_to_output_vis_image, ios::out | ios::binary | ios::trunc);
        if (output_vis_image.is_open())
        {
            // set pointer to the beginning of the image
            output_vis_image.seekg(0, ios::beg);
            for (unsigned int index_11 = 0; index_11 < no_pixels; index_11++) {
                output_vis_image.write((char*)&vis_image[index_11], sizeof(double));
            }
            output_vis_image.close();
        }
        else {
            std::cout << "Warning: Unable to open vis image file!!!" << "\n";
        }

        // stop to measure elapsed time at the end
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // print elapsed time in milliseconds, microseconds and nanoseconds
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[seconds]" << std::endl;
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[millisec]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[microsec]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[nanosec]" << std::endl;
    }

    // delete image buffer for foreground
    delete[] image_buffer_fg;
    // delete image buffer for background
    delete[] image_buffer_bg;
    // delete phase buffer for foreground
    delete[] phase_buffer_fg;
    // delete phase buffer for background
    delete[] phase_buffer_bg;
    // delete amplitude buffer for foreground
    delete[] amp_buffer_fg;
    // delete amplitude buffer for background
    delete[] amp_buffer_bg;
    // delete offset buffer for foreground
    delete[] offset_buffer_fg;
    // delete offset buffer for background
    delete[] offset_buffer_bg;
    // delete buffer for differential phase image
    delete[] dph_image;
    // delete buffer for absorption image
    delete[] abs_image;
    // delete buffer for visibility image
    delete[] vis_image;

    return 0;
}
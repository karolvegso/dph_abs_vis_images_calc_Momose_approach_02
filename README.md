# dph_abs_vis_images_calc_Momose_approach_02

This C++ simple program was written to perform calculation of differential phase images, absorption images nad visibility images for X-ray Talbot interferometry or grating interferometry in real time. The program opens M_fg = 5 foreground images and M_bg = 5 background images. The foreground and background images are encoded as 16 bit unsigned integers. The path to foreground folders is <path_to_fg_folder> and the path to background folders is <path_to_bg_folder>. The path to output folder is specified in <path_to_output_folder>. 

The constant in the program are number of steps in fringe scanning: M = 5, M_fg = 5, M_bg = 5. The size of input images is 1536 x 512 pixels. The foreground images are stored in 2D array <image_buffer_fg> and the background images are stored in the 2D array <image_buffer_bg>. 

The program was written for real time X-ray imaging. Therefore, it opens first subfolder with number <no_subfolder_fg_initial = 1> and the last subfolder with number <no_subfolder_fg_final = 11200> for foreground. The same is for background. It opens first subfolder with number <no_subfolder_bg_initial = 1> and the last subfolder with number <no_subfolder_bg_final = 11200> for background. The first subfolder for foreground and background is in the form 000001, and the last subfolder for foreground and background is in the form 011200. 

The result or final differential phase image is stored in 1D array of double values <dph_image>, the final absorbtion image is stored in the 1D array of double values <abs_image>, and the final visbility image is stored in the 1D array of double values <vis_image>. The differential phase, absorption and visbility images are saved to the output folder. 

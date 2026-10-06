Folder setup and Data
1.	The following folder contains data from 9 different single-volume acquisitions obtained from one eye. The data from these volumes are in folders 632, 634, 636, … , 650, respectively.
   
 <img width="468" height="150" alt="image" src="https://github.com/user-attachments/assets/e777c2f5-b1f7-4713-ad2d-dc72e26a8552" />
 
3.	Looking within each of these folders, we have each individual B-scan and the segmentation of the retina and choroid we used for each B-scan.
   
 <img width="468" height="202" alt="image" src="https://github.com/user-attachments/assets/3b15c87c-8e6f-45e2-a9ef-a2e6abe99640" />
 
At the end of each folder, there is a frame_3D.mat which contains a 3D matrix of the data, a projection image (en face), and non-OCT fundus images acquired during the acquisition.

<img width="468" height="139" alt="image" src="https://github.com/user-attachments/assets/d211593d-e7fb-4978-83b7-3ca3a7dc6837" />

 5.	As no standardized data type exists for different OCT devices, the user will need to determine how to acquire the structural B-scans and put the data into folders. Often, the OCT device manufacturer provides instructions for this.
6.	 We took a modified version of OCT Retinal Layer Segmenter to segment the choroid and retina. You can use any alternative segmentation method. The segmentation code we used is available here: https://github.com/FOIL-NU/RNFL_RPE_segmentation.

System Requirements
We used MATLAB 2023a to run the following code. 
We also borrowed the vesselness2D function for the Jerman filter used in point matching (original license for that function included).

Running the example
1.	Open Runme_example.m to run steps the steps of volumetric montaging after segmentation and data preparation. See comments for details.
2.	To adjust for different data, one needs to input the folder with the data.
   
<img width="438" height="100" alt="image" src="https://github.com/user-attachments/assets/7784fd0c-8c19-4265-af16-71e64c508710" />

4.	Next, the user needs to input which volumes are overlapping as well as assign names to files with the matching points.
   
 <img width="417" height="168" alt="image" src="https://github.com/user-attachments/assets/f33b7972-70a0-4c8b-a5ca-983a755b5ad9" />

6.	Finally, the user inputs which volume to use as the reference volume (usually 1, but the volume with the most overlapping volumes is recommended), the voxel size of the input volumes, and the voxel size of the output volumes.

<img width="417" height="171" alt="image" src="https://github.com/user-attachments/assets/eb6f5002-4627-4086-8482-3cd7bf881ee6" />

7.	The code will automatically identify matching points between the volumes.
   
 <img width="455" height="176" alt="image" src="https://github.com/user-attachments/assets/26e5dcad-c9c6-47d6-a56a-7fdc856987ba" />

The output points will be depicted on the en faces.

 <img width="468" height="244" alt="image" src="https://github.com/user-attachments/assets/c8d3a9a3-dd2d-40cb-b961-3f4e3f2ddd93" />

If the user dislikes the generation of matching points, they can run automatedmatching_quadrants_custompoints.m. See comments within these files for more details.
8.	Finally, we can run the function to volumetrically montage the volumes. F is a 3D matrix containing the information on the montaged volume. 

<img width="400" height="46" alt="image" src="https://github.com/user-attachments/assets/418fa624-ca1a-43b2-8685-97d8301a9167" />


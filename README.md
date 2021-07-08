# TissueSegmentationGUI
GUI for segmenting optical scans of ablated tissue

## Usage
The segmentation GUI can be used to create segmented masks from optical scans of ablated tissue. Scanned slices of tissue can be segmented automatically using a threshold, or manually via MATLAB's built in polygonal ROI creation tool. Tissue masks are then stored as 2D binary matricies where 1 indicates an ablated pixel, and 0 indicates an unablated pixel. These masks can then be collated and interpolated into a 3D mask using nearest neighbor, linear, spline or pchip interpolation. 
### Starting the GUI
On the main GUI page you will be presented with multiple inputs in which you will enter the geometry information of the 2D volumes as well as the parameters for the optical scans. Clicking 'Lock' will lock those parameters for the given trial and allow you to move onto the next stage of the GUI. 
### Adding data to the GUI
Press the select folder button and select a folder with the tissue scans you wish to segment. They will be loaded into memory and accessable via the drop down under the folder selection button. Navigating to a different scan can be done by selecting a file in this drop down. 
### Positioning
The GUI allows for translating and rotating the optical scan. The image can be translated via the two sliders on the side of the axes on the main window. Rotation is in the lowest panel of the GUI. Positive numbers will rotate clockwise, negative numbers will rotate counterclockwise. 
Once the tissue is aligned with the markers on the main axes change the slice number parameter to the desired slice number and press 'save image' to add it to the array of tissue images for further processing. 
### Segmentation
Once an image has been saved it can be segmented by going to the second tab group. Select the target slice via the dropdown on this page and click the 'draw rectangle' button to 

## Internals

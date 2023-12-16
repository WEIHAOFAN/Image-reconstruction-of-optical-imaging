function R1=StrongestContiguousRegion1(data,th,prior)

% This function segments out the strongest continuous region of non-zero
% values from a data set.  All other non-zero voxels are set to zero.
% Output volume has 1 contiguous region of activation with voxel values
% equal to input data values. 'th' is a threshold defined as a percentage
% of volume magnitude.
%
% If prior is defined, choose region which contains the prior.
%
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht & Jason Trobaugh
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


% Binarize data set
bw_data=data;
if exist('prior','var')
    MaxData=max(data(prior(1),prior(2),prior(3)));
    p=1;
else
    MaxData=max(data(:));
    p=0;
end
Threshold=MaxData*th;
bw_data=+(bw_data>=Threshold);


%% Separate largest contiguous region
L=bwlabeln(bw_data);     
stats=regionprops(L,'Area');
A=[stats.Area];

if p
    region=L(prior(1),prior(2),prior(3));
    bw_data(L~=region)=0;
else
    contribution=zeros(max(L(:)),1);
    for i=1:max(L(:))
        value=data(L==i);
        contribution(i)=sum(value(:));
    end
    [~,mI]=max(contribution);
    bw_data(L~=mI)=0;

end


%% Finalize output data
R1=data.*bw_data;
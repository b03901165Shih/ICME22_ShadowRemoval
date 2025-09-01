%   Input is a gray image of size ImgRow x ImgCol.
%   Output is a set of subimages, where each subimage is of size ImgRow x
%   ImgCol. These subimages are stacked up to form TextureFeatures (cell),
%   and the number of subimages = norient x nscale.
function TextureFeatures = OwnGetGaborResponse(GrayImg)

%   Predefined thresholds.
norient = 8;            %   Number of wavelet scales.
nscale = 5;             %   Number of filter orientations.
boundaryExtension = 32; %   Number of pixels to pad

orientationsPerScale = zeros(1,nscale);
orientationsPerScale(:) = norient;
ImgRow = size(GrayImg,1);   ImgCol = size(GrayImg,2);
GaborFilters = createGabor(orientationsPerScale, [ImgRow ImgCol] +2*boundaryExtension);

GrayImg = im2double(GrayImg);
% % scale intensities to be in the range [0 255]
% GrayImg = GrayImg - min(GrayImg(:));
% GrayImg = 255*GrayImg/max(GrayImg(:));
% % prefiltering: local contrast scaling
% %GrayImg = prefilt(GrayImg, 4);     %Alex Chia


TextureFeatures = OwnGaborConvolve(GrayImg, GaborFilters, boundaryExtension);
return;

end

% Input images are double in the range [0, 255];
% You can also input a block of images [ncols nrows 3 Nimages]
%
% For color images, normalization is done by dividing by the local
% luminance variance.
function output = prefilt(img, fc)

w = 5;
s1 = fc/sqrt(log(2));

% Pad images to reduce boundary artifacts
img = log(img+1);
img = padarray(img, [w w], 'symmetric');
[sn, sm, c, N] = size(img);
n = max([sn sm]);
n = n + mod(n,2);
img = padarray(img, [n-sn n-sm], 'symmetric','post');

% Filter
[fx, fy] = meshgrid(-n/2:n/2-1);
gf = fftshift(exp(-(fx.^2+fy.^2)/(s1^2)));
gf = repmat(gf, [1 1 c N]);

% Whitening
output = img - real(ifft2(fft2(img).*gf));
clear img

% Local contrast normalization
localstd = repmat(sqrt(abs(ifft2(fft2(mean(output,3).^2).*gf(:,:,1,:)))), [1 1 c 1]);
output = output./(.2+localstd);

% Crop output to have same size than the input
output = output(w+1:sn-w, w+1:sm-w,:,:);
return;

end



function TextureFeatures = OwnGaborConvolve(GrayImg, GaborFilters, boundaryExtension)

%   Doing filtering.

[ny nx Nfilters] = size(GaborFilters);
ImgRow = size(GrayImg,1);       ImgCol = size(GrayImg,2);
% pad image
GrayImg = padarray(GrayImg, [boundaryExtension boundaryExtension], 'symmetric');
GrayImg = single(fft2(GrayImg));


%   Convolving
TextureFeatures = zeros(ImgRow, ImgCol, Nfilters);
for n = 1:Nfilters
    ig = abs(ifft2(GrayImg.*repmat(GaborFilters(:,:,n), [1 1])));
    ig = ig(boundaryExtension+1:ny-boundaryExtension, boundaryExtension+1:nx-boundaryExtension, :);
    TextureFeatures(:,:,n)  = ig;
end

% %   Fast Convoling.
% TextureFeatures = zeros(ImgRow, ImgCol, Nfilters);
% TmpImg = zeros(size(GrayImg,1), size(GrayImg,2), Nfilters);
% for n = 1:Nfilters    
%     TmpImg(:,:,n) = GrayImg.*GaborFilters(:,:,n);
% end
% TmpImg = ifft2(TmpImg);
% TmpImg = abs(TmpImg);
% TextureFeatures = TmpImg(boundaryExtension+1:ny-boundaryExtension, boundaryExtension+1:nx-boundaryExtension, :);


return;

end

% Precomputes filter transfer functions. All computations are done on the
% Fourier domain. 
%
% If you call this function without output arguments it will show the
% tiling of the Fourier domain.
%
% Input
%     numberOfOrientationsPerScale = vector that contains the number of
%                                orientations at each scale (from HF to BF)
%     n = imagesize = [nrows ncols] 
%
% output
%     G = transfer functions for a jet of gabor filters
function G = createGabor(or, n)
Nscales = length(or);
Nfilters = sum(or);

if length(n) == 1
    n = [n(1) n(1)];
end

l=0;
for i=1:Nscales
    for j=1:or(i)
        l=l+1;
        param(l,:)=[.35 .3/(1.85^(i-1)) 16*or(i)^2/32^2 pi/(or(i))*(j-1)]; %#ok<AGROW>
    end
end

% Frequencies:
%[fx, fy] = meshgrid(-n/2:n/2-1);
[fx, fy] = meshgrid(-n(2)/2:n(2)/2-1, -n(1)/2:n(1)/2-1);
fr = fftshift(sqrt(fx.^2+fy.^2));
t = fftshift(angle(fx+sqrt(-1)*fy));

% Transfer functions:
G=zeros([n(1) n(2) Nfilters]);
for i=1:Nfilters
    tr=t+param(i,4); 
    tr=tr+2*pi*(tr<-pi)-2*pi*(tr>pi);

    G(:,:,i)=exp(-10*param(i,1)*(fr/n(2)/param(i,2)-1).^2-2*param(i,3)*pi*tr.^2);    
end


if nargout == 0
    figure
%     for i=1:Nfilters
%         max(max(G(:,:,i)));
%         contour(fftshift(G(:,:,i)),[1 .7 .6],'r');
%         hold on;
%         drawnow;
%     end
%     axis('on');
%     axis('square');
%     axis('ij');
    for i=1:Nfilters
        contour(fx, fy, fftshift(G(:,:,i)),[1 .7 .6],'r');
        hold on
    end
    axis('on')
    axis('equal')
    axis([-n(2)/2 n(2)/2 -n(1)/2 n(1)/2])
    axis('ij')
    xlabel('f_x (cycles per image)')
    ylabel('f_y (cycles per image)')
    grid on
    
    figure;
    for ii=1:Nfilters
        subplot(4,8,ii);
        w = fftshift(ifft2(G(:,:,ii)));
        w = real(w);  
        %w = (w - min(w(:))) / (max(w(:)) - min(w(:)));
        imshow(w/max(w(:)));
    end
    
    
for ii = 1:32
subplot(4,8,ii); 
end



end
return;

end



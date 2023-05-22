
function [imageHdr,hdrCompBitWidth] = HdrCombine(imageS, imageL, Color_BW, bitWidthImageL, times, dataOffset, transAreaStart, transAreaLength,maxBit,motionDetect)

if(Color_BW == 0)
    pattern = 'GRBG';
    imageL_G = DemosaicG(imageL,pattern,bitWidthImageL); % 将输入图像从RAW域转化到G空间
    imageS_G = DemosaicG((imageS.*times/64),pattern,bitWidthImageL); % 将输入图像从RAW域转化到G空间
else
    imageL_G = imageL;
end




widthExpand = 0;
transAreaLength = transAreaLength * 2^(bitWidthImageL-8);
transAreaEnd = min(256,transAreaStart) * 2^(bitWidthImageL-8) + transAreaLength;
transAreaStart = transAreaStart * 2^(bitWidthImageL-8);

transAreaWeight = zeros(1,2^bitWidthImageL); % 计算过渡区域的权重值
%transAreaWeight(1,transAreaEnd+1:2^bitWidthImageL) = 1;
transAreaWeight(1,transAreaEnd+1:2^bitWidthImageL) = 4096;
% kCoef = 1/(transAreaEnd-transAreaStart);
% kCoef = round(4096/(transAreaEnd-transAreaStart));

for i=transAreaStart:transAreaEnd
   kCoef = bitshift((i-transAreaStart),12);
   transAreaWeight(i) = fix(kCoef/(transAreaEnd-transAreaStart));
end

imageHdr = double(imageL);
imageS = max(0,double(imageS) - dataOffset);
%times = fix(expTimes * gainTimes * 64);

motionThSt = 64;
motionThLenExp = 6;
motionThEnd = 128;

MaxWeightMotion = 4096;


for i = 1:8
    if bitshift(times,-(i+6))==1
         widthExpand = i; 
         break;
    end
end
if(times-(bitshift(1,(widthExpand+6)))>0)
        widthExpand = widthExpand+1;
end
hdrCompBitWidth = ceil(bitWidthImageL + widthExpand);
imageLLog = log2(imageL)/12;
[rowNums,colNums] = size(imageS);
weight = zeros(rowNums,colNums);
weightL = zeros(rowNums,colNums);
weightMotion = zeros(rowNums,colNums);
weightNew = zeros(rowNums,colNums);
weightL3 = zeros(rowNums,colNums);
M = zeros(rowNums,colNums);

flag= zeros(rowNums,colNums);

WcLL = MertensContrast(imageL_G);
WcSS = MertensContrast(imageS_G);
Wc = zeros(rowNums,colNums);
Wc = (WcLL./(WcSS+WcLL+1))*4096;

se = strel('square',3);
B=[1 1 1 1 1;
   1 1 1 1 1;
   1 1 1 1 1;
   1 1 1 1 1;
   1 1 1 1 1];

B1 = [ 1 1 1;
   1 1 1;
   1 1 1];

% WcLLBin = WcLL>512;
%  WcLL = imerode(WcLL,B1);
WcLL = imdilate(WcLL,B);
  for i=2:rowNums-1
    for j=2:colNums-1       
        if(WcLL(i,j)  >= 1024)
           weightL3(i,j) = 0;
        else
            weightL3(i,j) = (1024-WcLL(i,j))*8;
        end
    end
  end


for i=2:rowNums-1
    for j=2:colNums-1
        if(i == 596) &&(j == 1045)
            aaaa= 1;
        end
        if(i == 537) &&(j == 496)
            aaaa= 1;
        end
        
          refVal = max(imageL_G(i,j),imageL(i,j));
         if(refVal>= transAreaStart && refVal< transAreaEnd)
            weightImageS = transAreaWeight(refVal);
            weightImageL = 4096 - weightImageS;
         elseif(refVal< transAreaStart)
            weightImageS = 0;
            weightImageL = 4096;
         else
            
            weightImageS = 4096;
            weightImageL = 0;
         end
%          weightImageLL = 4096-abs(bitshift(refVal,-7)-4096);
%          weightImageSS = 4096-abs(bitshift(imageS_G(i,j),-7)-4096);
% 
%          if ((weightImageLL*WcLL(i,j)+weightImageSS*WcSS(i,j)) == 0)
%              weightImageL = 0;
%          else 
%             weightImageL = floor(4096*weightImageLL*WcLL(i,j)/(weightImageLL*WcLL(i,j)+weightImageSS*WcSS(i,j)));
%          end
%             weightImageS = 4096-weightImageL;
%          
%          M(i,j) = ((weightImageLL*(refVal^2)+weightImageSS*(imageS_G(i,j)^2))*(weightImageLL+weightImageSS))/(weightImageLL*refVal+weightImageSS*imageS_G(i,j))^2;
%          M(i,j) = M(i,j)*4096-4096;
            
            
         diff = 0;
         for blki = -1:1
             for blkj = -1:1
                 pixS = floor((imageS(i+blki,j+blkj)*times)/64);
                 pixS = min(pixS,4096);
                 diff = diff + abs(imageL(i+blki,j+blkj)-pixS);
             end
         end
         
        diff = min(255,(diff/128));
        
        
% 
%         if(diff >= motionThEnd)
%            weightL2 = 0;
%         elseif(diff <= motionThSt)
%             weightL2 = MaxWeightMotion;
%         else
%             weightL2 = ((motionThEnd-diff)*MaxWeightMotion);
%             weightL2 = bitshift(weightL2,-motionThLenExp);
%         end
        

        if(diff >= motionThEnd)
           weightL2 = MaxWeightMotion;
        elseif(diff <= motionThSt)
            weightL2 = 0;
        else
            weightL2 = ((diff-motionThSt)*MaxWeightMotion);
            weightL2 = bitshift(weightL2,-motionThLenExp);
        end
        
        
%         if(WcLL(i,j)  == 50)
%            weightL3 = 0;
%         elseif(WcLL(i,j)  <= 10)
%             weightL3 = 4096;
%         else
%             weightL3 = ((50-WcLL(i,j) )*4096)/40;
% 
%         end
        weightL(i,j) = 4096-weightImageS;
        weightMotion(i,j) = weightL2;
%     if (maxBit == 20)         
%               weight(i,j) = min(min(weightL2,4096-weightImageS),weightL3(i,j));
% %              weight(i,j) = min(weightL2,4096-weightImageS);
%     else
%           weight(i,j) = min(weightL2,4096-weightImageS);
%     end
 
%weight(i,j) = min(weightL2,4096-weightImageS);

if (motionDetect)

       %weight(i,j) = min(weightL2,4096-weightImageS);
       weight(i,j) = 4096-((4096-weightL2)*weightImageS)/4096;
       if ((weightImageS ==4096)&&(weightL2==4096) )
        flag(i,j) =1;
       end
    
else
    weight(i,j) = 4096-weightImageS;
end

%         
%          weight(i,j) = weightImageL;
%          if (refVal>= transAreaStart)
%            
%           weight(i,j) = weightImageL;
%        else
%           weight(i,j) = 4096;
%        end


%          if ((4096-weightImageS)<4096)
%              
%             weightImageLL = 4096-abs(2*refVal-4096);
%             weightImageSS = 4096-abs(2*(imageS(i,j)*times)/64-4096);
%             if weightImageSS <0
%                 weightImageSS = 0;
%             end
% %             if (maxBit == 16)
% %                 weightImageSS = 4096-abs(bitshift((imageS(i,j)*times)/64,-3)-4096);
% %             elseif(maxBit == 20)
% %                 weightImageSS = 4096-abs(bitshift((imageS(i,j)*times)/64,-7)-4096);
% %             end
% 
%              if (weightImageLL + weightImageSS == 0)
%                  weightImageL = 0;
%              else 
%                 weightImageL = floor(4096*weightImageLL/(weightImageLL+weightImageSS));
%              end
%                 weightImageS = 4096-weightImageL;
%                 weight(i,j) = weightImageL;
%          else
%              weight(i,j) = 4096;
%          end
%          
%          if (refVal>transAreaEnd) &&(floor((imageS(i,j)*times)/64)>(transAreaEnd*16))
%              weight(i,j) = 2048;
%          end

    end   
end
 for i=1+4:rowNums-4
    for j=1+4:colNums-4
        BlkWeight = zeros(8,8);
        for a = -4:3
            for b = -4:3
                BlkWeight(a+5,b+5) = weight(i+a,j+b);
            end
        end
        denoise = floor(sum(sum(BlkWeight))/64);
        weightNew(i,j) = denoise;
        
    end
 end


 a = (2^maxBit - 2^(maxBit-4));
 
        

for i=3:rowNums-2
    for j=3:colNums-2

          % weight(i,j) = weightImageS;
            %weightImageL = 1 - weightImageS;
                 if (flag(i,j) == 1)
                    %imageL(i,j) = imageL(i,j)+((4096-weightL(i,j))*a/4096);
                    if((mod(i,2) == 1)&&(mod(j,2) == 1))%g
                        imageL(i,j) =4096;
                    elseif((mod(i,2) == 0)&&(mod(j,2) == 1))%B
                        imageL(i,j) =4096/2.4747;
                    elseif ((mod(i,2) == 1)&&(mod(j,2) == 0))%R
                        imageL(i,j) =4096/0.7536;
                    else
                        imageL(i,j) =4096;
                    end
                    
                   % imageL(i,j) =4096;
                end                
                weightL_imgL = floor(imageL(i,j) * weight(i,j) + 2048);
                
%                 if (flag(i,j) == 1)
%                     weightL_imgL = weightL_imgL*256;
%                 end
                
                
                weigthS_imgS = floor(imageS(i,j) * (4096-weight(i,j)) + 2048);
                imgS_times = bitshift(weigthS_imgS,-12) * times;
                imageHdr(i,j) = bitshift(weightL_imgL,-12) + bitshift(imgS_times,-6);
        
             imageHdr(i,j) =max(0,min(2^maxBit, imageHdr(i,j)));
            %imageHdr(i,j) = bitshift(imageHdr(i,j), -6);
            %imageHdr(i,j) = imageL(i,j) * weightImageL + imageS(i,j) * weightImageS * times;  
    end
end
            
end

function We = MertensWellExposedness(img, we_mean, we_sigma)
    %mean for the Well-exposedness weights.
    if(~exist('we_mean', 'var'))
        we_mean = 0.5; %as in the original  
    end
    %sigma for the Well-exposedness weights.
    if(~exist('we_sigma', 'var'))
        we_sigma = 0.2; %as in the original paper
    end
    sigma2 = 2.0 * we_sigma^2;
    [r, c, col] = size(img);
    We = ones(r, c);
    for i=1:col
        We = We .* exp(-(img(:,:,i) - we_mean).^2 / sigma2);
    end
end

function Wc = MertensContrast(L)
    H = [0 1 0; 1 -4 1; 0 1 0]; 
    Wc = abs(imfilter(L, H, 'replicate'));
end

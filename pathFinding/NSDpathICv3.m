function ICs = NSDpathICv3( SPorder, G, startPoint, ascFlag )

    for r=1:SPorder%
        ICs{r}(1)=startPoint; %easy one

        if SPorder>=2            %almost as easy one
            AscDescDir=1;
            if nargin == 4
                if ascFlag
                    AscDescDir=-1;
                end
            end

            %these are all the possible paths satisfying the SD DE
            ICs{r}(2)=(AscDescDir*1i*factorial(SPorder)/G{SPorder+1}(startPoint)).^(1/SPorder) * exp(1i*2*pi*(r)/SPorder);
        end

        for v=3:SPorder
            error('Have not coded this nasty looking FDB stuff yet');
        end
    end

end
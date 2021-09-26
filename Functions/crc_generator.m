function CRC = crc_generator(crcType)
   
    switch crcType
            case 'CRC24A'
                crcg =[1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1]; 
            case 'CRC24B'
                crcg = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];
            case 'CRC24C'
                crcg = [1 1 0 1 1 0 0 1 01 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
            case 'CRC16' 
                crcg = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
            case 'CRC11'
                crcg = [1 1 1 0 0 0 1 0 0 0 0 1];
            case 'CRC6'
                crcg = [1 1 0 0 0 0 1];
            case 'none'
                crcg = 1;
            otherwise
                error('%s is not a valid CRC type!',crcType);
    end
    
    CRC.g = crcg;
    CRC.L = length(crcg)-1;
end

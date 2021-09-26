
function crc = compute_crc(msg,crcg)
    crcL = length(crcg)-1;
    [~,rem] = gfdeconv([zeros(1,crcL) fliplr(msg)],fliplr(crcg));
    crc = fliplr([rem zeros(1,crcL-length(rem))]);
end
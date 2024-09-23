function code = TUDa_getColor(code)
    switch code
        case "0a"
            code = "#DCDCDC";
        case "1a"
            code = "#5D85C3";
        case "2a"
            code = "#009CDA";
        case "3a"
            code = "#50B695";
        case "4a"
            code = "#AFCC50";
        case "5a"
            code = "#DDDF48";
        case "6a"
            code = "#FFE05C";
        case "7a"
            code = "#F8BA3C";
        case "8a"
            code = "#EE7A34";
        case "9a"
            code = "#E9503E";
        case "10a"
            code = "#C9308E";
        case "11a"
            code = "#804597";
        case "0b"
            code = "#B5B5B5";
        case "1b"
            code = "#005AA9";
        case "2b"
            code = "#0083CC";
        case "3b"
            code = "#009D81";
        case "4b"
            code = "#99C000";
        case "5b"
            code = "#C9D400";
        case "6b"
            code = "#FDCA00";
        case "7b"
            code = "#F5A300";
        case "8b"
            code = "#EC6500";
        case "9b"
            code = "#E6001A";
        case "10b"
            code = "#A60084";
        case "11b"
            code = "#721085";
        case "0c"
            code = "#898989";
        case "1c"
            code = "#004E8A";
        case "2c"
            code = "#00689D";
        case "3c"
            code = "#008877";
        case "4c"
            code = "#7FAB16";
        case "5c"
            code = "#B1BD00";
        case "6c"
            code = "#D7AC00";
        case "7c"
            code = "#D28700";
        case "8c"
            code = "#CC4C03";
        case "9c"
            code = "#B90F22";
        case "10c"
            code = "#951169";
        case "11c"
            code = "#611C73";
        case "0d"
            code = "#535353";
        case "1d"
            code = "#243572";
        case "2d"
            code = "#004E73";
        case "3d"
            code = "#00715E";
        case "4d"
            code = "#6A8B22";
        case "5d"
            code = "#99A604";
        case "6d"
            code = "#AE8E00";
        case "7d"
            code = "#BE6F00";
        case "8d"
            code = "#A94913";
        case "9d"
            code = "#961C26";
        case "10d"
            code = "#732054";
        case "11d"
            code = "#4C226A";
        otherwise
            error("TUDa_get Color: color does not exist");
    end
end

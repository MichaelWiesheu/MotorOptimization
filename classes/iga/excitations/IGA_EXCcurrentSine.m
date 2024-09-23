classdef IGA_EXCcurrentSine < IGA_EXCcurrent
    properties (Access = public)
        IApp;
        Omega;
        UVWshift;
        Phase0;
    end

    methods
        function obj = IGA_EXCcurrentSine(parentElement, phasesStruct, iapp, omega, uvwshift, phase0, nWindings)
            currentFunction = @(t) iapp*sin(omega*t + uvwshift + phase0);
            obj@IGA_EXCcurrent(parentElement, phasesStruct, currentFunction, nWindings)
            obj.IApp = iapp;
            obj.Omega = omega;
            obj.UVWshift = uvwshift;
            obj.Phase0 = phase0;
            % Update current function with object values
            obj.Current = @(t) obj.IApp*sin(obj.Omega*t + obj.UVWshift + obj.Phase0);
        end
    end
end
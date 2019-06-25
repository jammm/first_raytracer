#!/bin/bash
set -eux

# xpra start :0 --no-pulseaudio

Xorg -noreset +extension GLX +extension RANDR +extension RENDER -logfile ./0.log -config /etc/X11/xorg.conf :0

/bin/bash


import datetime
import os
import threading

from cereal import log, messaging
from openpilot.common.params import Params
from openpilot.common.realtime import Priority, config_realtime_process
from openpilot.common.time import system_time_valid
from openpilot.system.hardware import HARDWARE

from openpilot.selfdrive.frogpilot.controls.frogpilot_planner import FrogPilotPlanner
from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_functions import FrogPilotFunctions, is_url_pingable
from openpilot.selfdrive.frogpilot.controls.lib.frogpilot_variables import FrogPilotVariables

WIFI = log.DeviceState.NetworkType.wifi

def automatic_update_check(started, params):
  update_available = params.get_bool("UpdaterFetchAvailable")
  update_ready = params.get_bool("UpdateAvailable")
  update_state_idle = params.get("UpdaterState", encoding='utf8') == "idle"

  if update_ready and not started:
    HARDWARE.reboot()
  elif update_available:
    os.system("pkill -SIGHUP -f system.updated.updated")
  elif update_state_idle:
    os.system("pkill -SIGUSR1 -f system.updated.updated")

def time_checks(automatic_updates, deviceState, now, started, params, params_memory):
  if not is_url_pingable("https://github.com"):
    return

  screen_off = deviceState.screenBrightnessPercent == 0
  wifi_connection = deviceState.networkType == WIFI

  if screen_off and wifi_connection:
    if automatic_updates:
      automatic_update_check(started, params)

def frogpilot_thread(frogpilot_toggles):
  config_realtime_process(5, Priority.CTRL_LOW)

  params = Params()
  params_memory = Params("/dev/shm/params")

  frogpilot_functions = FrogPilotFunctions()
  frogpilot_planner = FrogPilotPlanner()

  run_time_checks = False
  time_validated = system_time_valid()
  update_toggles = False

  pm = messaging.PubMaster(['frogpilotPlan'])
  sm = messaging.SubMaster(['carState', 'controlsState', 'deviceState', 'frogpilotCarControl', 'frogpilotCarState',
                            'frogpilotNavigation', 'frogpilotPlan', 'longitudinalPlan', 'modelV2', 'radarState'],
                            poll='modelV2', ignore_avg_freq=['radarState'])

  while True:
    sm.update()

    now = datetime.datetime.now()
    deviceState = sm['deviceState']
    started = deviceState.started

    if started and sm.updated['modelV2']:
      frogpilot_planner.update(sm['carState'], sm['controlsState'], sm['frogpilotCarControl'], sm['frogpilotCarState'],
                               sm['frogpilotNavigation'], sm['modelV2'], sm['radarState'], frogpilot_toggles)
      frogpilot_planner.publish(sm, pm, frogpilot_toggles)

    if FrogPilotVariables.toggles_updated:
      update_toggles = True
    elif update_toggles:
      FrogPilotVariables.update_frogpilot_params(started)

      update_toggles = False

    if now.second == 0:
      run_time_checks = True
    elif run_time_checks or not time_validated:
      time_check = threading.Thread(target=time_checks, args=(frogpilot_toggles.automatic_updates, deviceState, now, started, params, params_memory,))
      time_check.start()
      run_time_checks = False

      if not time_validated:
        time_validated = system_time_valid()
        if not time_validated:
          continue

def main():
  frogpilot_thread(FrogPilotVariables.toggles)

if __name__ == "__main__":
  main()

import logging
from logging import Logger


def create_logger(logger_config) -> Logger:
    logger = logging.getLogger(logger_config[LoggerKeys.NAME])

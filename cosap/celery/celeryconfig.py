import os

RABBITMQ_DEFAULT_USER = os.getenv("RABBITMQ_DEFAULT_USER")
RABBITMQ_DEFAULT_PASS = os.getenv("RABBITMQ_DEFAULT_PASS")

if not RABBITMQ_DEFAULT_USER or not RABBITMQ_DEFAULT_PASS:
    raise EnvironmentError(
        "RabbitMQ credentials not configured. "
        "Set RABBITMQ_DEFAULT_USER and RABBITMQ_DEFAULT_PASS environment variables."
    )

REDIS_URL = os.getenv("COSAP_REDIS_URL", "redis://redis:6379/0")

broker_url = (
    f"amqp://{RABBITMQ_DEFAULT_USER}:{RABBITMQ_DEFAULT_PASS}@rabbitmq:5672/vhost"
)
result_backend = REDIS_URL
task_track_started = True
worker_redirect_stdouts = False

# ---------- Task reliability defaults ----------
# Acknowledge tasks only after they complete (not when received).
# This ensures tasks are re-delivered if a worker crashes mid-execution.
task_acks_late = True

# If a worker is lost (killed, OOM, etc.), reject the task so the broker
# can re-deliver it to another worker.
task_reject_on_worker_lost = True

# Hard time limit: forcibly kill the worker process after this many seconds.
# Individual tasks can override via their decorator.
task_time_limit = int(os.getenv("COSAP_TASK_TIME_LIMIT", 172800))  # 48 hours

# Soft time limit: raise SoftTimeLimitExceeded so the task can clean up.
task_soft_time_limit = int(os.getenv("COSAP_TASK_SOFT_TIME_LIMIT", 169200))  # 47 hours

# Limit prefetch to 1 task per worker to avoid hogging long-running tasks.
worker_prefetch_multiplier = 1

# Serialization
task_serializer = "json"
result_serializer = "json"
accept_content = ["json"]
